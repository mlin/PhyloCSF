#!/bin/bash -e

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ] ; do SOURCE="$(readlink "$SOURCE")"; done
HERE="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

function test_in_project {
  dx-build-applet $HERE/.. -d $TEST_PROJECT:/
  dx="dx --project-context-id=$TEST_PROJECT"
  job=$($dx run :/PhyloCSF -i alignments="PhyloCSF development:/test data/tal-AA dm3_12flies alignments" \
          -y --brief)
  $dx wait $job
  $dx find jobs --origin $job
  $dx head $job:scores
}

if [ "$1" = "inner" ]; then
  test_in_project
  exit
fi

TEST_PROJECT=$(dx new project --brief "PhyloCSF test ($(date))")
export TEST_PROJECT

function cleanup {
  echo dx rmproject -y $TEST_PROJECT > /dev/null
}
trap "cleanup; exit" int

$SOURCE inner &
inner=$!

exit_code=0
wait $inner || exit_code=$?

if [ "$exit_code" -ne 0 ]; then
  echo "Leaving behind project for failed test, named \"$(dx describe $TEST_PROJECT --name)\""
  echo "To delete: dx rmproject -y $TEST_PROJECT"
else
  cleanup
fi

exit $exit_code
