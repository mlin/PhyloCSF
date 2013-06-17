#!/bin/bash -e

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ] ; do SOURCE="$(readlink "$SOURCE")"; done
HERE="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

function test_in_project {
  dx-build-applet $HERE/.. -d $TEST_PROJECT:/
  dx="dx --project-context-id=$TEST_PROJECT"
  job1=$($dx run :/PhyloCSF -i alignments="PhyloCSF development:/test data/tal-AA dm3_12flies alignments" \
           --delay-workspace-destruction -y --brief)
  job2=$($dx run :/PhyloCSF -i alignments="PhyloCSF development:/test data/tal-AA dm3_12flies alignments" \
           -i bls=true -i dna=true \
           --delay-workspace-destruction -y --brief)
  job3=$($dx run :/PhyloCSF -i alignments="PhyloCSF development:/test data/tal-AA dm3_12flies alignments" \
           -i frames=6 -i all_scores=true -i bls=true -i anc_comp=true -i dna=true -i aa=true \
           --delay-workspace-destruction -y --brief)
  job4=$($dx run :/PhyloCSF -i alignments="PhyloCSF development:/test data/ALDH2_exons hg19_29mammals alignments" \
           -i strategy=fixed -i orf=ATGStop -i frames=6 -i bls=true -i anc_comp=true -i aa=true \
           --delay-workspace-destruction -y --brief)
  job5=$($dx run :/PhyloCSF \
           -i alignments="PhyloCSF development:/test data/CCDS longest isoforms hg19_29mammals chr21 alignments" \
           -i strategy=fixed -i alignments_per_job=10 \
           --delay-workspace-destruction -y --brief)

  $dx wait $job1
  $dx export tsv -o - $job1:scores
  $dx wait $job2
  $dx export tsv -o - $job2:scores
  $dx wait $job3
  $dx export tsv -o - $job3:scores
  $dx wait $job4
  $dx export tsv -o - $job4:scores
  $dx wait $job5
  printf "%d/197\n" $($dx export tsv -o - --no-header $job5:scores | cut -f2 | grep -v - | wc -l)
}

if [ "$1" = "inner" ]; then
  test_in_project
  exit
fi

TEST_PROJECT=$(dx new project --brief "PhyloCSF test ($(date))")
export TEST_PROJECT

function cleanup {
  dx rmproject -y $TEST_PROJECT > /dev/null
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
