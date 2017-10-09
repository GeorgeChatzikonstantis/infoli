#!/bin/bash

HOME="$(getent passwd $USER | awk -F ':' '{print $6}')"
THREAD_FILE="$HOME/brainframe-reses/phaethon.state"
MICNAME="mic0"

ALL=`ssh $MICNAME "cat /proc/cpuinfo | grep "bogo" | wc -l"`
USED=`ssh $MICNAME "ps axHr| wc -l"`



echo $USED $ALL>$THREAD_FILE
echo $USED $ALL

