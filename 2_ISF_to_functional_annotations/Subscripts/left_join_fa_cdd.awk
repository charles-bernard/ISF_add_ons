#!/usr/bin/awk


BEGIN {
	FS="\t";
	file_idx = 0;
	k = 0;
}

FNR == 1 {
	file_idx = file_idx + 1;
}

####################################
# FIRST FILE: FASTA ################
####################################
file_idx == 1 && $0 ~ /^>/ {
	split($0, header_fields, /(>| )/);
	seq_id = header_fields[2];
	is[seq_id] = 0;
}

####################################
# SECOND FILE: OUTPUT OF RPS BLAST #
####################################
# Retrieve seq id and corresponding
# cdd id
file_idx == 2 && $0 !~ /^#/ {
	line[k] = $0;
	is[seq_id] = 1; 
	k++;
}

END {
	for(i = 0; i < k; i++) {
		print line[i];
	}

	for(s in is) {
		if(!is[s]) {
			fake_line = s;
			for(i = 0; i < NF; i++) {
				fake_line = fake_line "\t";
			}
			print fake_line
		}
	}
}

