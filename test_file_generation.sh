#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run run_aumic_search \
    python3 main.py ./config_file
assert_exit_code 0
assert_no_stderr
assert_in_stdout "Creating stellar FFDs"
assert_in_stdout "Starting on AU Mic"
assert_in_stdout "Operations for AU Mic finished."

declare -a dir_arr=(./data/ \
                    ./data/AU_Mic/ \
                    ./searches/)

declare -a file_arr=(./data/AU_Mic/AU_Mic_01.csv \
                     ./data/AU_Mic/AU_Mic_01_flares.ecsv \
                     ./data/AU_Mic/AU_Mic_FFD.png)

echo ""
echo "Checking if expected folders were generated..."
for dir in "${dir_arr[@]}"
do
    if [ -d "$dir" ]
    then
        echo "$dir exists."
    else
        echo "Could not find $dir."
    fi
done

echo ""
echo "Checking if expected files were generated..."
for file in "${file_arr[@]}"
do
    if test -f "$file"
    then
        echo "$file exists."
    else
        echo "Could not find $file."
    fi
done
