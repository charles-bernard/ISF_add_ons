#!/usr/bin/awk


BEGIN {
	FS="(\t)";
	file_idx = 0;
	toto = 0;

	one_letter[0] = "[X]"
	cog_supercat["[X]"] = "NO COG ASSOCIATED";
	cog_cat["[X]"] = "No COG associated";

 	k = 0;
 	l = 1;
}

FNR == 1 {
	file_idx = file_idx + 1;
}

####################################
# FIRST FILE: COG ONE LETTER CODE ##
####################################
# Retrieve the description associated
# with the letter
file_idx == 1 && $0 ~ /^[A-Z]/ {
	curr_super_category = $0;
}

file_idx == 1 && $0 ~ / \[[A-Z]\]/ {
	split($0, fields, " ");
	one_letter[l] = fields[1];
	description = fields[2];
	for(f = 3; f <= length(fields); f++) { 
		description = description " " fields[f];
	}
	cog_supercat[one_letter[l]] = curr_super_category;
	cog_cat[one_letter[l]] = description;
	l++;
}

####################################
# SECOND FILE: OUTPUT OF RPS BLAST #
####################################
# Retrieve seq id and corresponding
# cdd id
file_idx == 2 && $0 !~ /^#/ {
	seq_id[k] = $1;
	split($2, cdd_tag, "|");
	cdd_id[k] = cdd_tag[length(cdd_tag)];
	k++;
}

####################################
# THIRD FILE: CDD 2 COG CORRESP. ###
####################################
# retrieve correspondance between
# cdd_id and cog_id
file_idx == 3 && $2 ~ /^COG/ {
	cdd_cog[$1]=$2;
	cdd_description[$1] = $4;
}

####################################
# FOURTH FILE: COG ID 2 COG LETTER #
####################################
# retrieve correspondance between
# cog_id and cog_letter
file_idx == 4 && $0 ~ /^\[[A-Z]+\]/ {
	split($0, fields, " ");
	cog_id = fields[2];
	cog_letters = fields[1];
	cog_nb_letters[cog_id] = length(cog_letters) - 2 # because 2 brackets
	for(i = 0; i < cog_nb_letters[cog_id]; i++) {
		cog_letter[cog_id][i] = "[" substr(cog_letters, 2 + i, 1) "]"
	}
}

END {
	print "#Seq Id\tCOG Supercategory\tCOG Letter\tCOG Category\tCDD Id\tCDD description" \
		> comprehensive_out
	for(i = 0; i < k; i++) {
		curr_seq = seq_id[i];
		curr_cdd = cdd_id[i];

		if (!cdd_cog[curr_cdd]) {
			cdd_cog[curr_cdd] = "X";
		}
		curr_cdd_description = cdd_description[curr_cdd];

		curr_cog = cdd_cog[curr_cdd];
		curr_cog_nb_letters = cog_nb_letters[cog_id];

		for(j = 0; j < curr_cog_nb_letters; j++) {
			curr_cog_letter = cog_letter[curr_cog][j];
			curr_cog_cat = cog_cat[curr_cog_letter];
			curr_cog_supercat = cog_supercat[curr_cog_letter];

			print curr_seq "\t" curr_cog_supercat "\t" curr_cog_letter "\t" \
				curr_cog_cat "\t" curr_cdd "\t" curr_cdd_description \
				> comprehensive_out

			if(!occ[curr_cog_letter]) {
				occ[curr_cog_letter] = 1;
			} else {
				occ[curr_cog_letter] = occ[curr_cog_letter] + 1;
			}
		}
	}

	print "#COG Letter\tOccurence percentage in the family of genes\tCOG Supercategory\tCOG Category" \
		> summary_out;
	for(i = 0; i < l; i++) {
		curr_letter = one_letter[i];
		if(!occ[curr_letter]) {
			perc = 0;
		} else {
			perc = occ[curr_letter] / k * 100;
		}
		printf("%s\t%.2f%\t%s\t%s\n", curr_letter, perc, cog_supercat[curr_letter], cog_cat[curr_letter]) > summary_out;
	}
}