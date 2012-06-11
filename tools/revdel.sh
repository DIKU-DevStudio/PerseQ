#!/bin/sh
LOG=git.log
REV=revisions.log

i=0
for line in `cat $REV`;
do
    echo "Diff for $line (sprint $i)"
    git log --name-status $line > $LOG
    echo "Additions"
    grep '^A\s' $LOG | wc -l

    echo "Deletions"
    grep '^D\s' $LOG | wc -l

    echo "Revisions"
    grep '^M\s' $LOG | wc -l
    i=`expr $i + 1`
done

rm $LOG
