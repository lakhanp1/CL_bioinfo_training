use strict;
use warnings;
use LWP::Simple;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use URI::Escape;
use Pod::Usage;
use threads;
###############  Start Time  ########################################################
my $Time_Start="";
my $TS=time();
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";

##################  Main	 ######################################################

our $g_opts;
sub parse_opts{
	my $result = GetOptions(
					"feature_file|f=s" => \$g_opts->{'feature_file'},
					
					"start_column|s=i" => \$g_opts->{'start_column'},
					"end_column|e=i" => \$g_opts->{'end_column'},
					"chrom_column|c=i" => \$g_opts->{'chrom_column'},
					"direction_column|c=i" => \$g_opts->{'direction_column'},

					"start_range|sr=s" => \$g_opts->{'start_range'},
					"end_range|sr=s" => \$g_opts->{'end_range'},
					
					"start_bin|sb=i" => \$g_opts->{'start_bin'},
					"end_bin|eb=i" => \$g_opts->{'end_bin'},
					"start_step|ss=i" => \$g_opts->{'start_step'},
					"end_step|es=i" => \$g_opts->{'end_step'},
					
										
					"bin_count|bc=i" => \$g_opts->{'bin_count'},					
					"merge|m" => \$g_opts->{'merge'},
					"read_count|rc=i" => \$g_opts->{'read_count'},
					
					"socre_file|sf=s" => \$g_opts->{'socre_file'},
					
					"output_folder|of=s"=> \$g_opts->{'output_folder'},
					"outout_name|on=s"=> \$g_opts->{'outout_name'},
					
					"plot" => \$g_opts->{'plot'},
					"avg_plot" => \$g_opts->{'avg_plot'},
					"heatmap" => \$g_opts->{'heatmap'},
					
					#"threads=s"=> \$g_opts->{'threads'},
					


					"quiet"	 => sub { $g_opts->{'verbose'} = 0 },
					"help|?"	=> \$g_opts->{'help'}
				 );
	if(!($g_opts->{'feature_file'}) or !($g_opts->{'socre_file'}) ){
		print "ERR:01\n";
		&pod2usage( -verbose => 1);#exit status will be 1
	}
	
	if(!($g_opts->{'start_range'}) and !($g_opts->{'end_range'}) and !($g_opts->{'bin_count'})){
		print "ERR:02\n";
		&pod2usage( -verbose => 1);#exit status will be 1
	}
	if($g_opts->{'help'}){
		&pod2usage( -verbose => 1);#exit status will be 1
	}
}

&parse_opts();

my $feature_file=$g_opts->{'feature_file'};
my $socre_file=$g_opts->{'socre_file'};
#my $socre_type="counted";
#my $extension_length=-1;

my $chrom_column=0;
my $start_column=1;
my $end_column=2;

my $bin_count=20;
my $read_count=1;
$read_count=($g_opts->{'read_count'})/1000000 if(($g_opts->{'read_count'}));

my $start_range="-1000,0";
my $end_range="0,1000";

my $start_bin=50;
my $end_bin=50;
my $merge=0;

my $avg_plot="-";
my $heatmap="-";

my $start_r="";
my $end_r="";
my $bin_r="";

my $output_folder=dirname($feature_file);
my $outout_name=basename($feature_file)."-with-".basename($socre_file);

my $threads=1;


$start_column=$g_opts->{'start_column'}-1 if(($g_opts->{'start_column'}));

$end_column=$g_opts->{'end_column'}-1 if(($g_opts->{'end_column'}));

$chrom_column=$g_opts->{'chrom_column'}-1 if(($g_opts->{'chrom_column'}));

$start_range=$g_opts->{'start_range'} if(($g_opts->{'start_range'}));

$end_range=$g_opts->{'end_range'} if(($g_opts->{'end_range'}));

$start_bin=$g_opts->{'start_bin'} if(($g_opts->{'start_bin'}));

$end_bin=$g_opts->{'end_bin'} if(($g_opts->{'end_bin'}));

$bin_count=$g_opts->{'bin_count'} if(($g_opts->{'bin_count'}));

$output_folder=$g_opts->{'output_folder'} if(($g_opts->{'output_folder'}));

$outout_name=$g_opts->{'outout_name'} if(($g_opts->{'outout_name'}));

$start_r="-start-" if(($g_opts->{'start_range'}));

$end_r="-end-" if(($g_opts->{'end_range'}));

$bin_r="-rel-" if(($g_opts->{'bin_count'}));

$avg_plot="avg_plot" if(($g_opts->{'avg_plot'}));

$heatmap="heatmap" if(($g_opts->{'heatmap'}));

#$threads=$g_opts->{'threads'} if(($g_opts->{'threads'}));

my $start_step=$start_bin;
my $end_step=$end_bin;


print "$start_r\n";

$start_step=$g_opts->{'start_step'} if(($g_opts->{'start_step'}));
$end_step=$g_opts->{'end_step'} if(($g_opts->{'end_step'}));

my $direction_column=3;
$direction_column=$g_opts->{'direction_column'}-1 if(($g_opts->{'direction_column'}));

#$extension_length=$g_opts->{'extension_length'} if(($g_opts->{'extension_length'}));


print "\nParamater:\n\n";
print "feature_file: $feature_file\n";
print "socre_file: $socre_file\n";

print "chrom_column: $chrom_column\n";
print "start_column: $start_column\n";
print "end_column: $end_column\n";

print "start_range: $start_range\n" if(($g_opts->{'start_range'}));
print "end_range: $end_range\n" if(($g_opts->{'end_range'}));
print "start_bin: $start_bin bp\n" if(($g_opts->{'start_range'}));
print "end_bin: $end_bin bp\n" if(($g_opts->{'end_range'}));

print "bin_count: $bin_count\n" if(($g_opts->{'bin_count'}));
print "merge: True\n" if(($g_opts->{'merge'}));

print "read_count: $read_count\n" if(($g_opts->{'read_count'}));
print "output_folder: $output_folder\n";
print "outout_name: $outout_name\n";
#print "threads: $threads\n" if(($g_opts->{'threads'}));
print "R-Plot: True\n" if(($g_opts->{'plot'}));
print "avg_plot: True\n" if(($g_opts->{'avg_polt'}));
print "heatmap: True\n" if(($g_opts->{'heatmap'}));

print "\n";


my %score=();
my $kk=0;
open(IN,$socre_file)||die;
print "Reading $socre_file\n";
while(my $t=<IN>){
	chomp $t;
	next if($t=~/^track/);
	$t=~s/\r//;
	$t=~s/\n//;
	my @a=split(/\t/,$t);
	next if(!defined $a[1]);
	next if(!$a[1]=~/\d+/);
	if(defined $a[3]){
		$kk=1;
	}
	else{
		$kk=2;
	}
	last;
	
}
print "$kk\n";
if($kk==0){
	die "please check your score_file\n";
}

elsif($kk==1){
	open(IN,$socre_file)||die;
	print "Reading $socre_file\n";
	while(my $t=<IN>){
		chomp $t;
		next if($t=~/^track/);
		$t=~s/\r//;
		$t=~s/\n//;
		my @a=split(/\t/,$t);
		next if(!defined $a[1]);
		next if(!$a[1]=~/\d+/);
		for(my $i=$a[1];$i<=$a[2];$i++){
			$score{$a[0]}{$i}=$a[3];
		}
	
	}
}
else{
	my %score_tmp=();
	open(IN,$socre_file)||die;
	print "Reading $socre_file\n";
	my $s=0;
	my $ss=0;
	my $cc="";
	while(my $t=<IN>){
		chomp $t;
		next if($t=~/^track/);
		$t=~s/\r//;
		$t=~s/\n//;
		my @a=split(/\t/,$t);
		next if(!defined $a[1]);
		next if(!$a[1]=~/\d+/);
		if($cc ne $a[0]){
			$score{$cc}{$s}=$ss;
		}
		else{
			die "please check your score_file, it is should be the bed format or sorted sgr format\n" if($a[1]<$s);
			for(my $i=$s;$i<=$a[1];$i++){
				$score{$cc}{$i}=$ss;
			}
		}
		$s=$a[1];
		$cc=$a[0];
		$ss=$a[2];
	
	}
	
}
# open(IN,$socre_file)||die;
# print "Reading $socre_file\n";
# while(my $t=<IN>){
# 	chomp $t;
# 	my @a=split(/\t/,$t);
# 	next if(!defined $a[1]);
# 	next if(!$a[1]=~/\d+/);
# 	for(my $i=$a[1];$i<=$a[2];$i++){
# 		$score{$a[0]}{$i}=$a[3];
# 	}
# 	
# }


my @lines=();
open(IN,$feature_file)||die;
print "Reading $feature_file\n";
while(my $t=<IN>){
	chomp $t;
	$t=~s/\r//;
	$t=~s/\n//;
	next if($t=~/chr.*?\tstart\t/i);
	next if($t=~/\#/);
	next if($t=~/\!/);
	push(@lines,$t);
}


my $thread_arr=split_arr(\@lines,$threads);

my @threads_tmp=();
for(my $thr=0;$thr<$threads;$thr++){
	#my @tmp_arr=
	print "Thread:$thr\n";	
    $threads_tmp[$thr]=threads->create('WinSGR', \@{$$thread_arr[$thr]},\%score,$start_column,$end_column,$direction_column,$start_range,$end_range,$start_bin,$end_bin,$start_step,$end_step,$start_r,$end_r,$bin_r);
	
}

open(OUT1,">$output_folder/$outout_name.start.mat")||die;
open(OUT2,">$output_folder/$outout_name.end.mat")||die;
open(OUT3,">$output_folder/$outout_name.rel.mat")||die;
open(OUT4,">$output_folder/$outout_name.merge.mat")||die;

my %finial=();
for(my $thr=0;$thr<$threads;$thr++){
        my $res=$threads_tmp[$thr]->join();
    		foreach my $l (keys %{$$res{"start"}}){
    			$finial{"start"}{$l}=$$res{"start"}{$l};
    			#print $$res{"start"}{$l}."\n";
    		}

    		foreach my $l (keys %{$$res{"end"}}){
    			$finial{"end"}{$l}=$$res{"end"}{$l};
    		}
    		
    		foreach my $l (keys %{$$res{"rel"}}){
    			$finial{"rel"}{$l}=$$res{"rel"}{$l};
    			#print $$res{"rel"}{$l}."\n";
    		}
}

#die;

if($start_r ne ""){
	foreach my $l (@lines){
		my $str=$finial{"start"}{$l};
		print OUT1 "$l\t$str\n";
	}
}

if($end_r ne ""){
	foreach my $l (@lines){
		#my $str=join("\t",@{$finial{"end"}{$l}});
		my $str=$finial{"end"}{$l};
		print OUT2 "$l\t$str\n";
	}
}

if($bin_r ne ""){
	foreach my $l (@lines){
		#my $str=join("\t",@{$finial{"end"}{$l}});
		my $str=$finial{"rel"}{$l};
		print OUT3 "$l\t$str\n";
	}
}

if($end_r ne "" and $start_r ne "" and $start_r ne "" and ($g_opts->{'merge'})){
	foreach my $l (@lines){
		my $str1=$finial{"start"}{$l};
		
		my $str2=$finial{"rel"}{$l};
		
		my $str3=$finial{"end"}{$l};
		
		print OUT4 "$l\t$str1\t$str2\t$str3\n";
	}
}



if($start_r eq ""){
	system("rm $output_folder/$outout_name.start.mat");
}
if($end_r eq ""){
	system("rm $output_folder/$outout_name.end.mat");
}
if($bin_r eq ""){
	system("rm $output_folder/$outout_name.rel.mat");
}
if(!($g_opts->{'merge'})){
	system("rm $output_folder/$outout_name.merge.mat");
}

## R


# if(($g_opts->{'plot'})){
# 	if($start_r ne ""){
# 		system("rm $output_folder/$outout_name.start.mat");
# 	}
# 	if($end_r ne ""){
# 		system("rm $output_folder/$outout_name.end.mat");
# 	}
# 	if($rel_r ne ""){
# 		system("rm $output_folder/$outout_name.rel.mat");
# 	}
# 	if(($g_opts->{'merge'})){
# 		system("rm $output_folder/$outout_name.merge.mat");
# 	}
# }





sub WinSGR{
	my ($lines,$score,$start_column,$end_column,$direction_column,$start_range,$end_range,$start_bin,$end_bin,$start_step,$end_step,$start_r,$end_r,$bin_r)=@_;
	my $sub_arr=@$lines;
	#print "sub_arr:$sub_arr\n";
	#print "$start_column2,$end_column2,$start_range,$end_range,$start_bin,$end_bin,$start_step,$end_step,$start_r,$end_r\n";
	my %res=();
	if($start_r ne ""){
		print "start\n";
		my ($start_left,$start_right)=split(/,/,$start_range);
		#my $line_order=0;
		$start_step=int(abs($start_left-$start_right)/$start_bin);
		
		print "$start_step\n";
		
		foreach my $l (@$lines){
			#print "$l\n";
			my @a=split(/\t/,$l);
			my @res_arr=();
			if(defined $a[$direction_column]){
				if($a[$direction_column] eq "-"){
					for(my $i=$a[$end_column]-$start_left;$i>$a[$end_column]-$start_right;$i=$i-$start_bin){
						my $nn=0;
						my $sum=0;
						for(my $j=0;$j<$start_bin;$j++){
							last if($j+$i<$a[$end_column]-$start_right);
							$nn++;
							my $s=0;
							$s=$$score{$a[$chrom_column]}{($i+$j)} if(exists $$score{$a[$chrom_column]}{($i+$j)});
							$sum+=$s;
							#print "($i+$j)\t$$score{$a[0]}{($i+$j)}\t$sum\t$nn\n";
						}
						#print "$line_order\ti-$i\t\t\t$sum\t$nn\n";
						push(@res_arr,round($sum/($read_count*$nn),4)) if($nn ne "0");
						push(@res_arr,0) if($nn eq "0");
					}
				}
				else {
					for(my $i=$a[$start_column]+$start_left;$i<$a[$start_column]+$start_right;$i=$i+$start_bin){
						my $nn=0;
						my $sum=0;
						for(my $j=0;$j<$start_bin;$j++){
							last if($j+$i>$a[$start_column]+$start_right);
							$nn++;
							my $s=0;
							$s=$$score{$a[$chrom_column]}{($i+$j)} if(exists $$score{$a[$chrom_column]}{($i+$j)});
							$sum+=$s;
							#print "($i+$j)\t$$score{$a[0]}{($i+$j)}\t$sum\t$nn\n";
						}
						#print "$line_order\ti-$i\t\t\t$sum\t$nn\n";
						push(@res_arr,round($sum/($read_count*$nn),4)) if($nn ne "0");
						push(@res_arr,0) if($nn eq "0");
					}
				}
			}
			else{
				for(my $i=$a[$start_column]+$start_left;$i<$a[$start_column]+$start_right;$i=$i+$start_bin){
					my $nn=0;
					my $sum=0;
					for(my $j=0;$j<$start_bin;$j++){
						last if($j+$i>$a[$start_column]+$start_right);
						$nn++;
						my $s=0;
						$s=$$score{$a[$chrom_column]}{($i+$j)} if(exists $$score{$a[$chrom_column]}{($i+$j)});
						$sum+=$s;
						#print "($i+$j)\t$$score{$a[0]}{($i+$j)}\t$sum\t$nn\n";
					}
					#print "$line_order\ti-$i\t\t\t$sum\t$nn\n";
					push(@res_arr,round($sum/($read_count*$nn),4)) if($nn ne "0");
				push(@res_arr,0) if($nn eq "0");
				}
			}
			#$line_order++;
			#print "$a[$start_column]\t$start_left\t$start_right\n";
			
			
			#print "$line_order\t$l\t@res_arr\n";
			my $str=join("\t",@res_arr);
			#print "$str\n";
			$res{"start"}{$l}=$str;
		}
	}
	
	if($end_r ne ""){
		print "end\n";
		my ($end_left,$end_right)=split(/,/,$end_range);
		#my $line_order=0;
		$end_step=int(abs($end_left-$end_right)/$end_bin);
		foreach my $l (@$lines){
			my @a=split(/\t/,$l);
			my @res_arr=();
			if(defined $a[$direction_column]){
				if($a[$direction_column] eq "-"){
					#print "haha\t".($a[$start_column]+$end_left)."\t".($a[$start_column]-$end_right)."\t".$end_step."\n";
					for(my $i=$a[$start_column]+$end_left-1;$i>$a[$start_column]-$end_right;$i=$i-$end_bin){
						my $nn=0;
						my $sum=0;
						for(my $j=0;$j<$end_bin;$j++){
							last if($j+$i<$a[$start_column]-$end_right);
							$nn++;
							my $s=0;
							$s=$$score{$a[$chrom_column]}{($i+$j)} if(exists $$score{$a[$chrom_column]}{($i+$j)});
							$sum+=$s;
						}
						#print "$end_bin\n";
						push(@res_arr,round($sum/($read_count*$nn),4)) if($nn ne "0");
				push(@res_arr,0) if($nn eq "0");
					}
					
				}
				else {
					for(my $i=$a[$end_column]+$end_left+1;$i<=$a[$end_column]+$end_right;$i=$i+$end_bin){
						my $nn=0;
						my $sum=0;
						for(my $j=0;$j<$end_bin;$j++){
							last if($j+$i>=$a[$end_column]+$end_right);
							$nn++;
							my $s=0;
							$s=$$score{$a[$chrom_column]}{($i+$j)} if(exists $$score{$a[$chrom_column]}{($i+$j)});
							$sum+=$s;
						}
				
						push(@res_arr,round($sum/($read_count*$nn),4)) if($nn ne "0");
				push(@res_arr,0) if($nn eq "0");
					}
			
				}
			}
			else{
				for(my $i=$a[$end_column]+$end_left+1;$i<$a[$end_column]+$end_right;$i=$i+$end_bin){
					my $nn=0;
					my $sum=0;
					for(my $j=0;$j<$end_bin;$j++){
						last if($j+$i>$a[$end_column]+$end_right);
						$nn++;
						my $s=0;
						$s=$$score{$a[$chrom_column]}{($i+$j)} if(exists $$score{$a[$chrom_column]}{($i+$j)});
						$sum+=$s;
					}
				
					push(@res_arr,round($sum/($read_count*$nn),4)) if($nn ne "0");
				push(@res_arr,0) if($nn eq "0");
				}
			}
			
			#print "eee\t$line_order\t$l\t@res_arr\n";
			my $str=join("\t",@res_arr);
			$res{"end"}{$l}=$str;
		}
	}
	
	
	if($bin_r ne ""){
		print "in\n";
		#my $line_order=0;
		
		foreach my $l (@$lines){
			
			my @a=split(/\t/,$l);
			my ($rel_left,$rel_right)=($a[$start_column],$a[$end_column]);
			my $bin_size=int(abs($a[$start_column]-$a[$end_column])/$bin_count);
			#print "$bin_count\t$bin_size\t$rel_left,$rel_right,".($rel_right-$rel_left)."\n";

			#$line_order++;
			my @res_arr=();
			my $bin_tmp=0;
			
			for(my $i=0;$i<$bin_count;$i++){
				my $nn=0;
				my $sum=0;
				#print "j: ".($i*$bin_size)."\t".(($i+1)*$bin_size)."\t".($rel_left+$i*$bin_size)."\n";
				for(my $j=($i*$bin_size);$j<($i+1)*$bin_size;$j++){
					my $ind=$rel_left+$j;
					#print "$i\t$j\t$ind\t".($rel_left+$j)."\n";
					last if($ind>$a[$end_column]);
					$nn++;
					my $s=0;
					$s=$$score{$a[$chrom_column]}{(int($ind))} if(exists $$score{$a[$chrom_column]}{(int($ind))});
					#print "ok\n" if(exists $$score{$a[$chrom_column]}{(int($ind))});
					$sum+=$s;
				}
				#print "$sum\n";
				push(@res_arr,round($sum/($read_count*$nn),4)) if($nn ne "0");
				push(@res_arr,0) if($nn eq "0");

				
			}
			
			
			#print "eee\t$line_order\t$l\t@res_arr\n";
			if(defined defined $a[$direction_column]){
				if($a[$direction_column] eq "-"){
					@res_arr=reverse(@res_arr);
				}
			}
			my $str=join("\t",@res_arr);
			$res{"rel"}{$l}=$str;
		}
	}
	
	#print \%res_start."\t".\%res_end."\n";
	return \%res;
	
	
}


sub split_arr{
	my ($a,$threads)=@_;
	my $n=@$a;
	print "$n rows\n";
	my $sub_n=int($n/$threads)+1;
	my @res=();
	my $j=0;
	for(my $i=0;$i<$n;$i++){
		push(@{$res[$j]},$$a[$i]);
		$j++ if(($i+1)%$sub_n==0);
	}
	return \@res;
}

sub round{
	my ($s,$c)=@_;#$c, how many decimals want to keep
	my $n1=10**$c;
	#print "$c\t$n1\n";
	return $s=int(($s*$n1)+0.5)/$n1;
	
}


###############  End Time  ########################################################
my $Time_End="";
my $TE=time();
my $UseTime=$TE-$TS;
my $danwei="Seconds";
my $UT="";
if($UseTime>120 and $UseTime<7200){
	my $m=int($UseTime/60);
	my $s=$UseTime%60;
	$UT= "Used Time: $m minutes $s seconds";
}
elsif($UseTime<=120){
	$UT= "Used Time: $UseTime seconds";
}
elsif($UseTime>=7200){
		my $h=int($UseTime/3600);
		my $m=int(($UseTime%3600)/60);
		my $s= ($UseTime%3600)%60 ;
		$UT= "Used Time: $h hours $m minutes $s seconds";
}
 $Time_End = sub_format_datetime(localtime(time()));
 print "\nEnd Time :[$Time_End]\nTime-consuming: $UT\n";

###############################################################################
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
		$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


__END__

=head1 NAME

script - Using Getopt::Long and Pod::Usage

=head1 SYNOPSIS

script [options] [args ...]

Options: 

	-feature_file|f            * Feature_file or chromosome coordinate file, format:"chr  start  end"
	-socre_file|sf             * socre_file, like .sgr .bed(counted)

	-chrom_column|s            Int. In feature_file, which column is chromosome name, counting from 1, Default is 1
	-start_column|s            Int. In feature_file, which column is start site, counting from 1, Default is 2
	-end_column|e              Int. In feature_file, which column is end site, counting from 1, Default is 3
	-direction_column|e        Int. In feature_file, which column is direction, counting from 1, Default is 4
	
	-start_range|sr            Defined which range be scanned, Default is -1000,0
	-end_range|sr              Defined which range be scanned, Default is 0,1000
	
	-start_bin|sb              Int. Window size for scanning start_range, units is bp, Default is 50 
	-end_bin|eb                Int. Window size for scanning end_range, units is bp, Default is 50   
	-start_step|ss             Int. Step size for scanning start_range, units is bp, Default is 50
	-end_step|es               Int. Step size for scanning end_range, units is bp, Default is 50
	
	#for relative range
	-bin_count|b               Int. How many bins you want to check, Default is 20   
	-merge                     Bool. whether merge start, end, and rel score together, Default is False
	-read_count|rc             Int. read_count for normalization the score , Default is not to do normalization.
	
	#for output
	-output_folder|of          Output result folder, Default is same folder with feature_file
	-outout_name|on            Output result name, Default is same feature_file name combine socre_file name
	

	
	-help			           brief help message
	
	                          * mean option required

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut


Regards
Zhengqiang Miao

zqmiao@umac.mo
miaozhengqiang1987@gmail.com


