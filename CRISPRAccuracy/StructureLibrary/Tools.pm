package StructureLibrary::Tools;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(read_file_return_array
			write_file_from_array
			convert_array_to_string
			recursively_return_files_in_dir
			round_Xdp
			get_path_base_suffix
			);
@EXPORT_OK = qw(plot_results_with_R
			);


########################################################################################
# This method takes a file name and returns each line in an array, all spaces
# on the end of the line are removed.
# input: name of the file to be read.
# output: the array of all lines in the file.
########################################################################################
sub read_file_return_array{
	my ($file) = @_;
	
	open(IN_HANDLE, "<$file") || die "couldn't open the following file in package Tool,".
									 " sub read_file: $file/n";
		chomp(my @rows = <IN_HANDLE>);
	close(IN_HANDLE);
	
	return @rows;
}


##########################################################################################
# This method takes an array, which contains a string for each element that represents
# one line in the output file and writes the information in the array to  a file with the 
# given name.
# Input: 
#	filename 		The name of the file
#	output_array	Contents to write to the file, line by line
#	append			Append to existing file
##########################################################################################
sub write_file_from_array{
	my($filename, $output_array, $append) = @_;
	
	if($append){
		open(OUT_HANDLE, ">>$filename") || die "couldn't open the file for writing in package Tools,".
								"sub write_file_from_array: $filename\n";
	} else{
	
		open(OUT_HANDLE, ">$filename") || die "couldn't open the file for writing in package Tools,".
								"sub write_file_from_array: $filename\n";
	}
		
	chomp(@{$output_array});						
	foreach my $line (@{$output_array}){
		
		print OUT_HANDLE $line."\n";
	}
	
#	$, = "\t";
#	$\ = "\n";
#	print LIST_HANDLE @{$featurelist_ref};
#	$, = "";
#	$\ = "";
	
	close(OUT_HANDLE);
	return 1;
}

###########################################################################################
# This method creates a box plot in a pdf by reading a file with the results
# Input: null
# Output: null
#########################################################################################
sub plot_results_with_R {
	
	my $file = shift;
		
	# define an output file name for the boxplot in pdf format
	my $plot = "intron_exon_length_statistic.pdf";
	
	# give out information of what the script is doing
	print "plotting box plots and save it in: $plot\n";
	
	# opening filehandle and piping in the R command, which means everthing written into the
	# filehandle is carried out in R. 
	open(RCMD, "| R --vanilla --slave");
	
	# load the library called ggplot, this has to be manually installed first
	print RCMD "library(\"ggplot2\")\n";
	
	# read in file
	print RCMD "mydata<-read.table(\"$file\", header=T)\n";
	
	# create the boxplots and save it in "ieplot"
	print RCMD "myboxplot <- qplot(class, variable, data=mydata, geom=\"boxplot\",", 
			"xlab=\"Class\", ylab=\"Variable measure\", ",
			"main=\"Title of boxplot\")\n";
	
	# save the plot in the designated pdf
	print RCMD "ggsave(\"$plot\",plot=ieplot)\n";
	
	# close filehandle
	close(RCMD);
	
	# remove the files you don't need afterwards, cleanup.
	# it is possible to carry out system commands using system("some command")
	# you have to be aware of the fact that the commands differ between operating systems
	system("rm Rplots.pdf");
	return 1;
}

###################################################################################
# converts the given array into one-lined string with the given separator between
# the elements.
# Input: array, separator symbol (optional, otherwise default used)
# Output: the array as a string
##################################################################################
sub convert_array_to_string{
	my ($aref, $sep) = @_;
	
	if (!$sep){
		$sep = "\t";
		print "WARNING: no separator defined, so setting it to default, \\t\n"	
	}
	
	my $string = "";
	my $size = @{$aref};
	for (my $i = 0 ; $i < $size-1 ; $i++) {
		$string .= $aref->[$i].$sep;
	}
	$string .= $aref->[-1];
	return $string;
	return 1;
}


###################################################################################
# Returns all files in directory and all subdirectories. You can also specify
# a file prefix or part of the file name, so that only these are returned.
# 
# Input: 
#	path		The directory path
#	file_suffix	The suffix of files to be returned
#	file_subname The substring of the name of the files to be returned
# Output: 
#	An array of file names with their paths.
##################################################################################
sub recursively_return_files_in_dir {
    my ($path, $file_suffix, $file_subname) = @_;
    
    my $FUNCTION = "recursively_return_files_in_dir in package Tools";

    # Open the directory.
    opendir (DIR, $path)
        or die "ERROR in $FUNCTION: Unable to open $path: $!";

    # Read in the files.
    my @files = 
    	    	# add the path to the filename
    			map { $path . '/' . $_ }
    			# remove files '.' and '..'
    			grep { !/^\.{1,2}$/ } 
    			# read all files from the directory into an array
				readdir (DIR);
    

    # Close the directory.
    closedir (DIR);
    
    my @selected_files;
    
	# go through each file in given directory
    for (@files) {

        # If the file is a directory
        if (-d $_) {
            # Call the function again
            my @return_files = recursively_return_files_in_dir ($_, $file_suffix, $file_subname);
            push(@selected_files, @return_files);

		# add file to selected files
        } else { 
        	# add the files with the given suffix
        	if($file_suffix && ($_ =~ /\.?$file_suffix$/)){
        		push (@selected_files, $_);
        	# add the files that include the given subname
        	} elsif ($file_subname && ($_ =~ /$file_subname/)){
        		push (@selected_files, $_);
        	# add all files
        	} else {
        		push (@selected_files, $_) unless ($file_suffix || $file_subname);
        	}
        }
    }
    
    # return selected files
    return @selected_files;
}

##################################################################################
# Round the given number to x decimal places.
# Input
#	num 	The number to round
#	dp 		The number of decimal places
# Output 
# The rounded number
##################################################################################
sub round_Xdp{
	my ($num, $dp) = @_;
	my $potenz = 1;
	for(my $i=0;$i<$dp;$i++){
		$potenz *= 10;
	}
	my $n = $num*$potenz;
	if($n-int($n)< 0.4999999999999){
		return (int ($n)) / $potenz;
	} else {
		return (int ($n+1)) / $potenz;
	}
}

##################################################################################
# This method returns an array of information about the given file.
# The first element contains the path, the second element contains the
# file base name and the third element is the file suffix.
# Input
#	file	The file name with complete path
# Output
#	(path, base, suffix) An array with the file path in the first position,
#	the base of the file name until the last dot, and finally the file suffix
#	after the final dot.
##################################################################################
sub get_path_base_suffix{
	my $file = shift;
	
	die "There is no file name given!\n" unless ($file);
	my $file_len = length($file);
	
	my ($path, $base, $suffix) = "";
	my $finalslash = rindex($file,"/");
	my $finaldot = rindex($file,".");
	
	$path = substr($file, 0, $finalslash);
	$base = substr($file, ($finalslash+1), ($finaldot - $finalslash - 1));
	$suffix = substr($file, ($finaldot+1), ($file_len - $finaldot - 1));
	
	return ($path, $base, $suffix);
}

##################################################################################
# Calculates the median of an array of numbers
# Input
#	numbers 	The array of numbers
# Output 
# The median number of the input array
##################################################################################
sub median {
	my (@numbers) = @_;
	
	#Prevent division by 0 error in case you get junk data
	return undef unless(scalar(@numbers));
	
	my @sorted_numbers =  sort {$a <=> $b} @numbers;
	my $total = @numbers;
	my $mid = $total/2;
	my $pos1 = int($mid);
	my $pos2 = -1;
	## even
	if($pos1*2 == $total){
		$pos2 = $val1 + 1;
		return (($sorted_numbers[$pos1 - 1] + $sorted_numbers[$pos2 - 1])/2);
		
	## odd
	} else {
		return ($sorted_numbers[$pos1 - 1]);
	}
	return -1;
}

##################################################################################
# Calculates the standard deviation of an array of numbers
# Input
#	numbers 	The array of numbers 
# Output 
# An array with first the standard deviation and then the mean, so this won't need
# to be re-calculated if it is needed as well.
##################################################################################
sub standard_deviation {
	my(@numbers) = @_;
	#Prevent division by 0 error in case you get junk data
	return undef unless(scalar(@numbers));
	
	# Step 1, find the mean of the numbers
	my $total1 = 0;
	foreach my $num (@numbers) {
		$total1 += $num;
	}
	my $mean1 = $total1 / (scalar @numbers);
	
	# Step 2, find the mean of the squares of the differences
	# between each number and the mean
	my $total2 = 0;
	foreach my $num (@numbers) {
		$total2 += ($mean1-$num)**2;
	}
	my $mean2 = $total2 / (scalar @numbers);
	
	# Step 3, standard deviation is the square root of the
	# above mean
	my $std_dev = sqrt($mean2);
	return ($std_dev, $mean1);
}



#### HOWTO ####
# Sorting numerically: sort {$a<=>$b} @array