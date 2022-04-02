#!/bin/bash

qstat -t -u naim |sed '1,5d' |awk '$10=="R"' >current_jobs.txt
while [ -s current_jobs.txt ]; do
  while IFS= read -r line; do
    t1=$(echo "$line" |awk '{print $9}')
    t2=$(echo "$line" |awk '{print $11}')
    jobnum=$(echo "$line" |awk '{print $1}')
    t1hour=$(echo $t1 |cut -d":" -f1)
    t1min=$(echo $t1 |cut -d":" -f2 |sed 's/^0*//g')
    t2hour=$(echo $t2 |cut -d":" -f1)
    t2min=$(echo $t2 |cut -d":" -f2 |sed 's/^0*//g')
    if [ $t1hour -eq $t2hour ]; then
      diff=`expr $t1min - $t2min`
      if [ $diff -lt 2 ]; then
         echo "Killing ${jobnum}"
         qdel $jobnum
      fi
    fi
  done < current_jobs.txt
  sleep 5
  qstat -t -u naim |sed '1,5d' |awk '$10=="R"' >current_jobs.txt
done
echo "No jobs found"

