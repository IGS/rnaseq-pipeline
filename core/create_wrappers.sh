#!/usr/bin/bash

# create_wrappers.sh - Create wrapper files based on the executable scripts present in the bin directory
# Args - [1] - Base path of Ergatis package

script_dir=$(dirname "$0")

perl_exe=$(which perl)

for file in `ls -1 $1/bin`; do
	fileext=${file##*.}
	case $fileext in
		pl)
			$perl_exe ${script_dir}/perl2wrapper_ergatis.pl INSTALL_BASE=$1 $file
			;;
		py)
			$perl_exe ${script_dir}/python2wrapper_ergatis.pl INSTALL_BASE=$1 $file
			;;
		jl)
			$perl_exe ${script_dir}/julia2wrapper_ergatis.pl INSTALL_BASE=$1 $file
			;;
		*)
			echo "$file is not a Perl, Python, or Julia file... skipping\n"
			;;
	esac
done

# give 555 permissions to wrappers and executables
chmod 555 $1/bin/*
