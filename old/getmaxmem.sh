memory=`ps auxww |grep java |awk '{print $6}' |sed '1d'`
max=2000000
processes=`ps auxww |grep java |wc -l`
while [ processes -gt 1 ]
do
  if [ max -lt memory ]
  then
    memory=max
  fi
  memory=`ps auxww |grep java |awk '{print $6}' |sed '1d'`
done
