#!/usr/bin/env perl
###
# This script generated as described at http://uswest.ensembl.org/info/data/biomart/biomart_perl_api.html
###
# Some dependencies:
#  apt install libdbi-perl
#  apt install liblog-log4perl-perl
#  apt install libxml-dom-perl
###
# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict;
no warnings 'uninitialized';
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my $confFile = "/home/app/lib/perl/biomart-perl/conf/martURLLocation.xml";

##  For Biomart Central Registry navigate to http://www.biomart.org/biomart/martservice?type=registry
#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

my $action='clean';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

$query->setDataset("hsapiens_gene_ensembl");
$query->addAttribute("ensembl_gene_id");
$query->addAttribute("ensembl_gene_id_version");
$query->addAttribute("ensembl_transcript_id");
$query->addAttribute("ensembl_transcript_id_version");
$query->addAttribute("ensembl_peptide_id");
$query->addAttribute("ensembl_peptide_id_version");
$query->addAttribute("external_gene_name");

$query->formatter("TSV");

my $query_runner = BioMart::QueryRunner->new();

############################## GET COUNT ############################
#$query->count(1);
#$query_runner->execute($query);
#print $query_runner->getCount();
#####################################################################


############################## GET RESULTS ##########################
# to obtain unique rows only
$query_runner->uniqueRowsOnly(1);

$query_runner->execute($query);
$query_runner->printHeader();
$query_runner->printResults();
$query_runner->printFooter();
#####################################################################
