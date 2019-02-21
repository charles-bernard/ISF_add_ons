#!/usr/bin/awk


BEGIN {
	FS="\t";
	file_idx = 0;
	k = 0;

	# thresh variable is sent to this script
}

FNR == 1 {
	file_idx = file_idx + 1;
}

####################################
# FIRST FILE: SEED IDS #############
####################################
file_idx == 1 {
	is_seed[$1] = 1;
}

####################################
# SECOND FILE: EDGES IT1 ###########
####################################
# Retrieve seq id and corresponding
# cdd id
file_idx == 2 {
	query_species = $1;
	pident = $3;
	target_species = $11;

	if(!is_seed[target_species] && pident > thresh) {
		print $0;
	}
}


