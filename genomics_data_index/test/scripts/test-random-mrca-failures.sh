#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_PROJECT_DIR=$DIR/../../../

if [ "$#" -ne 1 ];
then
	echo "Usage: $0 [ATTEMPTS]"
	echo "Synopsys: This is a script to attempt to narrow-down random failures related to cluster scoring"
	echo "          specifically, realted to calculating mrca. This runs the Python test code in a loop"
	echo "          [ATTEMPTS] times and counts the number of failures."
	exit 0
fi

attempts=$1
attempt=0
success=0
fail=0

cd $ROOT_PROJECT_DIR
echo "Running for $attempts attempts"
while [ "$attempt" -lt "$attempts" ];
do
	echo -n "Attempt $attempt ... "
	pytest $ROOT_PROJECT_DIR/genomics_data_index/test/integration/api/query/test_ClusterScorer.py::test_score_samples_with_mutation_tree 2>/dev/null 1>/dev/null
	ret_value=$?
	if [ "$ret_value" = 0 ];
	then
		success=$((success+1))
		echo "SUCCESS"
	else
		fail=$((fail+1))
		echo "FAIL"
	fi

	attempt=$((attempt+1))
done
cd $DIR

percent=`echo "scale=1; 100 * ${fail}/${attempts}" | bc`
echo -e "Failed: ${percent}% (${fail}/${attempts})"
