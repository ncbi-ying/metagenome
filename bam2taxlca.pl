#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw"$Bin";
my $usage=<<USAGE;
bam2taxlca.pl : statistics for taxonomy and lca form database alignment bam
Usage:
	perl $0 <bam> <outprefix>
USAGE
die $usage unless @ARGV == 2;
my ($bam,$pref) = @ARGV;
my $samtools = "software_dir/samtools";

my $taxformat = "DB.tax.format";
#GenomeID	TaxName	TaxID	Kingdom	Type	Species	SpeTaxID	SpeGroup	SpeGroupTaxID	Genus
#GCA_000010525.1	Azorhizobium caulinodans	7	Bacteria	G-	Azorhizobium caulinodans	7	-	-	Azorhizobium	6	5369772

my $plist = "all.sort.info";
#GenBank_ID	GenBank_info	GenomeID	TaxID	Genome_Length
#AP009384.1	AP009384.1 Azorhizobium caulinodans ORS 571 DNA, complete genome	GCA_000010525.1	7	5369772

my $process = 4;
my %seq2gid; # seqid -> genomeid
my %tax2king;	# taxonomy -> kingdom
my %tax2tid;	# taxonomy -> taxid
my %tax2type;	# taxonomy -> type
my %taxtree;	# taxonomy tree
my %tinfo;	# taxonomy information full
open INF,"<",$plist or die $!;
while (<INF>) { # load information for seqid -> gnomeid
	chomp;
	my @a = split /\t/;
	$seq2gid{$a[0]} = $a[2];
}
close INF;
open TX,"<",$taxformat or die $!;
while (<TX>) { # load information for taxonomy dastabase
	chomp;
	my @a = split /\t/; #GenomeID       TaxName TaxID   Kingdom Species SpeTaxID        SpeGroup        SpeGroupTaxID   Genus   GenTaxID
	for my $idx (qw"5 7 9") {
		next if exists $tax2tid{$a[$idx]};
		my $idx1 = $idx + 1;
		$tax2tid{$a[$idx]} = $a[$idx1];
		$tax2king{$a[$idx]} = $a[3];
		$tax2type{$a[$idx]} = $a[4];
	}
	$tinfo{$a[0]} = \@a;
}
close TX;
#open IN,"$samtools view -@ $process $bam |" or die $!;
`$samtools sort -n -@ $process $bam -o $pref.sort.bam`;
open IN,"$samtools view $pref.sort.bam |" or die $!;
open my $taxfh,">","$pref.tax" or die $!;
open my $lcafh,">","$pref.tax.lca" or die $!;
print $taxfh "#ReadID\tTaxID\tKingdom\tGenus\tComplex\tSpecies\tMU\n";
print $lcafh "#ReadID\tTaxID\tKingdom\tGenus\tComplex\tSpecies\n";
my $lastid = '';
my %genu;	# genus uniq count
my %comu;	# complex uniq count
my %speu;	# species uniq count
my %spea;	# species aligned all
my %tax;	# taxonomy tempory storage
$tax{'gen'} = {};
$tax{'com'} = {};
$tax{'spe'} = {};
while (<IN>) {
	chomp;
	my @ll = split /\t/;
	my $sid = $ll[0];
	my $rid = $ll[2];
	next if $rid eq "*"; #20200414 gjp
	my $gid = $seq2gid{$rid};
	my @info = @{$tinfo{$gid}};
	next if ($info[5] eq '-' and $info[7] eq '-' and $info[9] eq '-');
	$lastid ||= $sid;
	if ($lastid ne $sid) { # summary for last readid
		&stat_readalign(\%tax,\%genu,\%comu,\%speu,\%spea,\%tax2king,\%tax2tid,$taxfh,$lcafh,$lastid);
		$lastid = $sid;
		$tax{'gen'} = {};
		$tax{'com'} = {};
		$tax{'spe'} = {};
	}
	$tax{'gen'}{$info[9]}++;
	$tax{'com'}{$info[7]}++;
	$tax{'spe'}{$info[5]} .= join("\t",$sid,$info[6],$info[3],$info[9],$info[7],$info[5])."\n";
	$taxtree{$info[3]}{$info[9]}{$info[7]}{$info[5]}++;
	$genu{$info[9]}+=0;
	$comu{$info[7]}+=0;
	$speu{$info[5]}+=0;
	$spea{$info[5]}+=0;
}

&stat_readalign(\%tax,\%genu,\%comu,\%speu,\%spea,\%tax2king,\%tax2tid,$taxfh,$lcafh,$lastid);

open TAB,">","$pref.anno" or die $!;
print TAB "#TaxID\tKingdom\tGram+/-\tGenus\tGenCount\tComplex\tComCount\tSpecies\tSpeUniq\tSpeMulti\tMU\n";
foreach my $king (qw"Bacteria Fungi Viruses Parasite") { # not output Archaea
	foreach my $genid ( keys %{$taxtree{$king}} ) {
		$genu{$genid} ||= 0; # intitialize genus for sort
	}
	foreach my $genid ( sort {$genu{$b}<=>$genu{$a}} keys %{$taxtree{$king}} ) {
		my $genc = $genu{$genid} || 0;
		foreach my $comid ( keys %{$taxtree{$king}{$genid}} ) {
			$comu{$comid} ||= 0; # initialize complex for sort
		}
		my %spe_u;
		my %spe_a;
		foreach my $comid ( keys %{$taxtree{$king}{$genid}} ) {
			if ($comid ne '-' and $comu{$comid}) {
				$spe_u{"$comid:-"} = $comu{$comid};
				$spe_a{"$comid:-"} = $comu{$comid};
			}
			foreach my $speid ( keys %{$taxtree{$king}{$genid}{$comid}} ) {
				next unless ($spea{$speid});
				$spe_u{"$comid:$speid"} = $speu{$speid} || 0;
				$spe_a{"$comid:$speid"} = $spea{$speid};
			}
		}
		foreach my $outid ( sort { $spe_u{$b}<=>$spe_u{$a} or $spe_a{$b}<=>$spe_a{$a} } keys %spe_u ) {
			my ($comid,$speid) = (split /:/,$outid);
			if ($speid eq '-') {
				my $ctax = $tax2tid{$comid};
				my $type = $tax2type{$comid};
				print TAB "$ctax\t$king\t$type\t$genid\t$genc\t$comid\t$comu{$comid}\t-\t-\t-\t-\n";
			} else {
				my $ctax = $tax2tid{$speid};
				my $suniq = $speu{$speid} || 0;
				my $sall = $spea{$speid} || 0;
				my $spem = $sall - $suniq;
				if ($speid eq '-') {
					$suniq = '-';
					$spem = '-';
				}
				$genc ||= '-';
				my $comc = $comu{$comid} || '-';
				my $type = $tax2type{$speid};
				my $mu = 0;
				if ($suniq == 0) {
					$mu = 9999;
				} else {
					$mu = $spem/$suniq;
				}
				print TAB "$ctax\t$king\t$type\t$genid\t$genc\t$comid\t$comc\t$speid\t$suniq\t$spem\t$mu\n";
			}
		}
	}
}
close TAB;
exit 0;

sub stat_readalign {
	my ($tax,$genu,$comu,$speu,$spea,$tax2king,$tax2tid,$taxfh,$lcafh,$lastid) = @_;
	my $gn = (keys %{$tax{'gen'}});
	my $genus = '-';
	if ($gn == 1 and (keys %{$tax{'gen'}})[0] ne '-') { # only one genus
		$$genu{(keys %{$tax{'gen'}})[0]}++;
		$genus = (keys %{$tax{'gen'}})[0];
	}
	my $complex = '-';
	my $cn = (keys %{$tax{'com'}});
	if ($cn == 1 and (keys %{$tax{'com'}})[0] ne '-') { # only one complex
		$$comu{(keys %{$tax{'com'}})[0]}++;
		$complex = (keys %{$tax{'com'}})[0];
	}
	my $species = '-';
	my $sn = (keys %{$tax{'spe'}});
	my $mu = 'multi';
	if ($sn == 1 and (keys %{$tax{'spe'}})[0] ne '-') { # only one species
		$$speu{(keys %{$tax{'spe'}})[0]}++;
		$mu = 'uniq';
		$species = (keys %{$tax{'spe'}})[0];
	}
	foreach my $spe (keys %{$tax{'spe'}}) { # set all statistics
		next if ($spe eq '-');
		$$spea{$spe}++;
		my @rinfo = split /\n/,$$tax{'spe'}{$spe};
		foreach my $cinfo (@rinfo) {
			print $taxfh "$cinfo\t$mu\n";
		}
	}
	if ($genus ne '-' or $complex ne '-' or $species ne '-') { # skip all '-' annotations for lca
		my $kingdom;
		my $ctaxid;
		if ($species ne '-') {
			$kingdom = $$tax2king{$species};
			$ctaxid = $$tax2tid{$species};
		} elsif ($complex ne '-') {
			$kingdom = $$tax2king{$complex};
			$ctaxid = $$tax2tid{$complex};
		} elsif ($genus ne '-') {
			$kingdom = $$tax2king{$genus};
			$ctaxid = $$tax2tid{$genus};
		}
		print $lcafh "$lastid\t$ctaxid\t$kingdom\t$genus\t$complex\t$species\n";
	}
	return 0;
}
