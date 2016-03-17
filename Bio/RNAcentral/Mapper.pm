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

    Utilities for mapping low-level INSDC accessions to the toplevel genomic
    coordinates using Ensembl API.

=cut


package Bio::RNAcentral::Mapper;

use strict;
use warnings;

use Exporter 'import';

our @EXPORT_OK = qw(&get_toplevel_mapping);
our @EXPORT = qw(&get_toplevel_mapping);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);


=head2 get_toplevel_mapping

    Arg [1]    : Reference to Ensembl Slice Adaptor
    Arg [2]    : INSDC accession
    Arg [3]    : Starting coordinates relative to the INSDC accessions
    Arg [4]    : Ending coordinates relative to the INSDC accessions
    Arg [5]    : Sequence
    Description: Map low-level accession to toplevel genomic coordinates
    Return     : seq_region_name (chromosome, patch etc), start, end, strand.

=cut

sub get_toplevel_mapping {

    my ($slice_adaptor, $accession, $local_start, $local_end, $sequence) = @_;

    my ($seq_region_name, $start, $end, $strand);

    # swap start and end if necessary
    if ($local_start > $local_end) {
        my $temp = $local_start;
        $start = $local_end;
        $end = $temp;
    }

    if ( my $slice = $$slice_adaptor->fetch_by_region( 'toplevel', $accession ) ) {

        my $feat = new Bio::EnsEMBL::Feature(
              -slice  => $slice,
              -start  => $local_start,
              -end    => $local_end,
              -strand => 1, # always use 1, but then confirm directionality using the actual sequence
        );

        if ( my $toplevel = $feat->transform('toplevel') ) {
            # try using transform to get toplevel coordinates
            $start = $toplevel->start();
            $end = $toplevel->end();
            $seq_region_name = $toplevel->seq_region_name();
            $strand = get_strand($sequence, $toplevel);
        } else {
            # if transform didn't work, try using project
            # this can return slices on different coordinate systems
            # e.g. chromosome and patch coordinates
            my $projection = $feat->project('toplevel');
            foreach my $segment ( @{$projection} ) {
                my $to_slice = $segment->to_Slice();

                $start = $to_slice->start();
                $end = $to_slice->end();
                $seq_region_name = $to_slice->seq_region_name();
                $strand = get_strand($sequence, $to_slice);

                if ($local_end - $local_start + 1 != $to_slice->length) {
                    print "Incorrect length\n";
                }
            }

            # my $display_id = $slice->display_id();
            # print "Errored on $display_id\n";
        }

    }

    # discard results if the strand couldn't be determined
    if (defined($strand) && $strand == 0) {
        $seq_region_name = $start = $end = $strand = undef;
    }

    return ($seq_region_name, $start, $end, $strand);
}


=head2 get_strand

    Arg [1]    : RNAcentral sequence, DNA form
    Arg [2]    : Ensembl Slice object
    Description: Given a genomic slice and a sequence,
                 find whether the sequence comes from
                 the forward or the reverse strand.
    Return     : strand. Possible values:
                 +1
                 -1
                 0 if sequence was not found in the slice.

=cut

sub get_strand() {

    my ($sequence, $slice) = @_;
    my $strand;

    my $rev = $sequence;
    reverse_comp(\$rev);

    if (index($sequence,$slice->seq) != -1 && index($rev,$slice->seq) != -1) {
        # sequence found on both strands
        $strand = 0;
    } elsif (index($sequence,$slice->seq) != -1) {
        $strand = $slice->strand();
    } elsif (index($rev,$slice->seq) != -1) {
        $strand = - $slice->strand(); # reverse complement matched the sequence, so reverse the strand of the feature
    } else {
        print "sequence not found\n";
        $strand = 0;
    }

    return $strand;
}


1;
