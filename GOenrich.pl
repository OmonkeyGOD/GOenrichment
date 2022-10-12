#use lib '/home/wangy660/perl5/lib';
use Math::CDF;
use Array::Utils qw(:all);
my %hash_GO; #GO -> genes
my %hash_gene; #gene -> GOs
my @arr_gene;
my %hash_GO_anno; #GO -> anno
my %hash_GO_ratio; #P(C) inferred based on reference set there are 5907 genes in L.elogisporus and 3849 have GO annotations. We should use 5907
my $tmp = "LELG_03853";
my @homo_high;
my $xx;

open X, "homo_high";
for $xx(<X>){
	chomp($xx);
	push(@homo_high, $xx);
	
}
close(X);
open F, "Lodderomyces_anno.txt";
open O, ">>res.txt";
for $line(<F>){
	chomp($line);
    @rec=split("\t", $line);
    if($rec[0] eq $tmp){
        push(@{$hash_gene{$rec[0]}}, $rec[1]);
    }else{
        push(@arr_gene,$tmp);
        #print $tmp,"\n";
        $tmp = $rec[0];
        push(@{$hash_gene{$rec[0]}}, $rec[1]);
    }
    if (!exists $hash_GO_anno{$rec[1]}) {
        $hash_GO_anno{$rec[1]}=$rec[2];
    }
    push(@{$hash_GO{$rec[1]}}, $rec[0]);
}
push(@arr_gene,$tmp);




foreach $key(keys(%hash_gene)){
    
    #print $key, "=>",@{$hash_gene{$key}},"\n";
    #print $key, "=>",$hash_gene{$key}[0],"\n";
}
my $gene_number = 3849; ##!#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@# 
foreach $key(keys(%hash_GO)){
	$num = @{$hash_GO{$key}};
	$hash_GO_ratio{$key} = $num / $gene_number; # P(C)
	#print $hash_GO_ratio{$key} ,"\n";
	#print O $key, "=>",@{$hash_GO{$key}},"\n";
    	#print $key, "=>",$hash_gene{$key}[0],"\n";
}


$len=@arr_gene;
foreach $i(0..$len-1){
    
    #print $arr_gene[$i],"\n";
}
close(F);

my @arr_gff;
#my $sca="NW_001813677.1";
print "sca please";
my $sca = <STDIN>;
chomp $sca;
my @arr_sca; #all genes in target sca
open G, "GCF_000149685.1_ASM14968v1_genomic.gff";
for $d(<G>){
    chomp($d);
    @arr_gff=split("\t", $d);
    if($arr_gff[0] eq $sca and $arr_gff[2] eq "gene"){
        $arr_gff[8]=~/ID=gene-(.*);Dbxref/;
        push(@arr_sca, $1);
    }
}
close(G);
my $flag;
print O $gene_number,';', $sca, "\n";
print O "GO\tGOlabel\tpvalue\t+/-represented\n";
foreach $key(keys(%hash_GO_anno)){
	my @arr_sca_updated = intersect(@arr_gene,@arr_sca); # only consider the genes with annotations
	#my @arr_sca_updated = @arr_sca; #consider all genes in scaffold no matter whether they have annotations
	#my @arr_sca_updated = intersect(@homo_high,@arr_gene);
	#my @arr_sca_updated = @homo_high;
	my @isect = intersect(@{$hash_GO{$key}},@arr_sca_updated);
	#if (@isect != 0){
		
		my $n = @arr_sca_updated;
		my $C = @isect;
		#$num = @{$hash_GO{$key}}; #tmp
		#$num_e=$num - $C; ##tmp
		#$hash_GO_ratio{$key} = $num_e / ($gene_number-$n); ##tmp
        	my $pvalue = &Math::CDF::pbinom($C,$n,$hash_GO_ratio{$key});
	
	#print "ok","\n";
	if($pvalue <= 0.05){
		if (($hash_GO_ratio{$key} * $n) > $C){
			$flag = '-';
			#print $flag;
		}else{
			$flag = '+';
		}
		$genes=join(",",@isect);
		$aaaa=@{$hash_GO{$key}};
		print  "number_in_list,list_size,GOsize","\n";
		print  $C , "\t",$n,"\t",$aaaa,"\t",$hash_GO_ratio{$key},"\n";
		print  $key,"\t",$hash_GO_anno{$key},"\t",$pvalue,"\t",$flag,"\t",$genes,"\n";
	#}
        #print @{$hash_GO{$key}},"\n";
    }
    #print $key, "=>",$hash_GO_anno{$key},"\n";
    #print $key, "=>",$hash_gene{$key}[0],"\n";
}

close(O);
