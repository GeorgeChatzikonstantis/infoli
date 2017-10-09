#!/bin/bash

HOME="$(getent passwd $USER | awk -F ':' '{print $6}')"
THREAD_FILE="$HOME/brainframe-reses/phaethon.state"
MICNAME="mic0"

ALL=`ssh $MICNAME "cat /proc/cpuinfo | grep "bogo" | wc -l"`
USED=`ssh $MICNAME "ps axHr| wc -l"`



echo $ALL $USED>$THREAD_FILE
echo $ALL $USED

read VAR1 VAR2 < $THREAD_FILE

((VAR1=$VAR1 + $VAR2))
echo $VAR1
echo $VAR2
