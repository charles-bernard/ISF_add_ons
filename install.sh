#!/bin/bash

ISF_ADDONS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";

# Change owner
# chown -R $USER:$USER $ISF_ADDONS_DIR

# checking Python version
if [[ $(which python3) != "/usr/bin/python3" ]]; then
    echo "No python3 found: please contact your system administrator";
    exit;
fi

#####################################
# INSTALLATIONs #####################
#####################################
function update {
	update_status="$1";
	if [[ $update_status==0 ]]; then
		apt-get update;
	fi
	update_status=1;
}

update_status=0;
mkdir "$ISF_ADDONS_DIR"/dependencies
cd "$ISF_ADDONS_DIR"/dependencies

# install blast2 if necessary
if ! [ -x "$(command -v rpsblast)" ]; then
	update $update_status;
	apt-get install blast2;
fi

# install gnu-parallel if necessary
parallel_version=0;
if [ -x "$(command -v parallel)" ]; then
	parallel_version=`parallel --version | head -n 1 | egrep -o '20[0-9]{2}'`;
fi
if (( $parallel_version < 2018 )); then
	wget -q http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2;
	tar -xvjf parallel-latest.tar.bz2;
	rm parallel-latest.tar.bz2;
	cd parallel*/;
	./configure;
	make;
	make install;
	cd ../;
fi

# install cdd2cog
# wget https://raw.githubusercontent.com/aleimba/bac-genomics-scripts/master/cdd2cog/cdd2cog.pl

cd ../

#####################################
# DOWNLOADs #####################
#####################################
mkdir "$ISF_ADDONS_DIR"/required_files
cd "$ISF_ADDONS_DIR"/required_files

# download cddid
wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
gunzip cddid.tbl.gz

# download cog files
wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt
wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog

cd ../
