#!/bin/bash

GLOBAL_LOG=jobs_cur.log
APPEND_LOG=jobs_all.log

echo "data_in,reference,jobnumber,owner,hostname,qsub_time,start_time,end_time,ru_wallclock,ru_utime,ru_stime,cpu,mem,maxvmem,failed,exit_status,command" > $GLOBAL_LOG

if [ ! -f $APPEND_LOG ]; then
	echo "data_in,reference,jobnumber,owner,hostname,qsub_time,start_time,end_time,ru_wallclock,ru_utime,ru_stime,cpu,mem,maxvmem,failed,exit_status,command" > $APPEND_LOG
fi

for file in `ls *.o*`
do
	echo Parsing filename $file...
	JID=`echo $file | awk 'BEGIN {FS = "\\\.o"}{print $2}'`
	echo Processing job $JID...
	LOG=$JID.log

	if [ -f $JID.cmd ]; then
		CMD=`cat $JID.cmd`
	else
		CMD="unknown"
	fi

	if [ -f $JID.data_in ]; then
		DATAIN=`cat $JID.data_in`
	else
		DATAIN="unknown"
	fi

	if [ -f $JID.ref_in ]; then
		REFIN=`cat $JID.ref_in`
	else
		REFIN="unknown"
	fi

	qacct -j $JID > $LOG

	HOST=`awk 'NR==3 {print $2;exit}' $LOG`
	OWNER=`awk 'NR==5 {print $2;exit}' $LOG`
	SUB_TIME=`awk 'NR==13 {print $2" "$3" "$4" "$5" "$6;exit}' $LOG`
	START_TIME=`awk 'NR==14 {print $2" "$3" "$4" "$5" "$6;exit}' $LOG`
	END_TIME=`awk 'NR==15 {print $2" "$3" "$4" "$5" "$6;exit}' $LOG`
	FAILED=`awk 'NR==18 {print $2;exit}' $LOG`
	EXIT_STATUS=`awk 'NR==19 {print $2;exit}' $LOG`
	TIME_CLOCK=`awk 'NR==20 {print $2;exit}' $LOG`
	TIME_USER=`awk 'NR==21 {print $2;exit}' $LOG`
	TIME_SYS=`awk 'NR==22 {print $2;exit}' $LOG`
	CPU_USE=`awk 'NR==38 {print $2;exit}' $LOG`
	MEM_USE=`awk 'NR==39 {print $2;exit}' $LOG`
	MAX_MEM=`awk 'NR==42 {print $2;exit}' $LOG`

	echo $DATAIN,$REFIN,$JID,$OWNER,$HOST,$SUB_TIME,$START_TIME,$END_TIME,$TIME_CLOCK,$TIME_USER,$TIME_SYS,$CPU_USE,$MEM_USE,$MAX_MEM,$FAILED,$EXIT_STATUS,$CMD >> $GLOBAL_LOG
	echo $DATAIN,$REFIN,$JID,$OWNER,$HOST,$SUB_TIME,$START_TIME,$END_TIME,$TIME_CLOCK,$TIME_USER,$TIME_SYS,$CPU_USE,$MEM_USE,$MAX_MEM,$FAILED,$EXIT_STATUS,$CMD >> $APPEND_LOG
done
