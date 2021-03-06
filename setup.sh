#!/bin/bash
# Copyright [2009-2014] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


########################################################
## Setup working environment for the Ensembl Perl API ##
########################################################

# get the directory of this script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# set the environment for Ensembl
PERL5LIB=${PERL5LIB}:${DIR}/bioperl-live
PERL5LIB=${PERL5LIB}:${DIR}/ensembl/modules
export PERL5LIB

# load environment variables with sensitive data
source ${DIR}/scripts/params.sh
