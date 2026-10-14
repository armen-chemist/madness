#!/bin/bash
# Stop ALL queues started with run_queue.sh and all running magneticatom programs of this user.
pkill -u "$USER" -f "run_queue.sh" && echo "queues stopped" || echo "no queue running"
sleep 1
pkill -u "$USER" -x magneticatom && echo "magneticatom stopped" || echo "no magneticatom running"
