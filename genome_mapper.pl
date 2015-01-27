#!/usr/bin/env perl

=pod

=head1 LICENSE

    Copyright [2009-2014] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

=head1 DESCRIPTION

    Script for mapping low-level INSDC accessions to the toplevel genomic
    coordinates using Ensembl API.

=cut


use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::RNAcentral::Mapper qw(get_toplevel_mapping);

use DBI;
use DBD::Oracle;
use Getopt::Long;


# initialisation
# the environmental variables are set in setup.sh
my $species;
my $help = '';

if ( !GetOptions( 'species|s=s' => \$species,
                  'help|h!'     => \$help )
     || !( defined $ENV{'ORACLE_USER'} &&
           defined $ENV{'ORACLE_PASSWORD'} &&
           defined $ENV{'ORACLE_SID'} &&
           defined $ENV{'ORACLE_PORT'} &&
           defined $ENV{'ORACLE_HOST'} )
     || $help
     || ! defined $species )
{
  print <<END_USAGE;
Usage:
  $0 --species=species
  $0 --help
    --species / -s  Species scientific name.
    --help    / -h  To see this text.
Example usage:
  $0 -s Homo_sapiens
END_USAGE

  exit(1);
}

# replace underscores in species names
$species =~ tr/_/ /;

# connect to the RNAcentral database
my $dsn = sprintf("dbi:Oracle:host=%s;service_name=%s;port=%s",
                  $ENV{'ORACLE_HOST'}, $ENV{'ORACLE_SID'},
                  $ENV{'ORACLE_PORT'});

my $dbh = DBI->connect($dsn, $ENV{'ORACLE_USER'}, $ENV{'ORACLE_PASSWORD'});

# enable Oracle error reporting
$dbh->func(1000000, 'dbms_output_enable');

# enable reading large CLOBs
$dbh->{LongReadLen} = 66000;
$dbh->{LongTruncOk} = 1;

# connect to Ensembl API
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_multiple_dbs(
    {-host => 'ensembldb.ensembl.org', # Ensembl
     -port => 5306,
     -user    => 'anonymous'
    },
    {-host => 'mysql-eg-publicsql.ebi.ac.uk', # Ensembl Genomes
     -port => 4157,
     -user => 'anonymous'
    }
);

my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );

# get locations to be mapped from the RNAcentral database
my $cmd = <<SQL;
    SELECT t1.accession as accession, primary_accession, local_start, local_end, strand,
           t3.seq_short, t3.seq_long, t2."NAME" as CHROMOSOME, primary_start, primary_end
    FROM rnc_accessions t1, rnc_coordinates t2, rna t3, xref t4
    WHERE
        t1.accession = t2.accession AND
        t4.ac = t2.accession AND
        t3.upi = t4.upi AND
        t1.species = '$species'
SQL

my $query_sql = $dbh->prepare($cmd);
$query_sql->execute();

# map low-level INSDC accessions to toplevel genomic coordinates
my $found = 0;
my $not_found = 0;

# prepare an SQL statement for saving the mappings
$cmd = <<SQL;
    UPDATE rnc_coordinates
    SET primary_start = ?, primary_end = ?, name = ?, strand = ?
    WHERE
        accession = ? AND
        Local_start = ? AND
        local_end = ?
SQL

my $update_sql = $dbh->prepare_cached($cmd);

my $i=0;
# iterate over the INSDC accessions, map them to their genomic coordinates,
# and save the results in the database
while (my $data = $query_sql->fetchrow_hashref()) {

    my ($seq_region_name, $start, $end, $seq, $strand);

    # # temporary limit
    # if ($i > 10) {
    #     $query_sql->finish();
    #     last;
    # } else {
    #     $i++;
    # }

    $seq = $data->{'SEQ_SHORT'} || $data->{'SEQ_LONG'};

    printf(
        "%s %i %i %i\n",
        $data->{'PRIMARY_ACCESSION'}, $data->{'LOCAL_START'},
        $data->{'LOCAL_END'}, $data->{'STRAND'}
    );

    ($seq_region_name, $start, $end, $strand) = get_toplevel_mapping(
         \$slice_adaptor,
         $data->{'PRIMARY_ACCESSION'},
         $data->{'LOCAL_START'},
         $data->{'LOCAL_END'},
         $seq);

    if ($seq_region_name) {
        $found++;
        $update_sql->execute($start, $end, $seq_region_name, $strand,
                             $data->{'ACCESSION'}, $data->{'LOCAL_START'},
                             $data->{'LOCAL_END'}) or warn('Error!');
    } else {
        $not_found++;
    }
}

# clean up
print "Mapped $found, not mapped $not_found\n";

$dbh->disconnect();
