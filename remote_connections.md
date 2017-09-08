### Tips and tricks for connecting to Roberts Lab computers remotely.


#### Tmux to keep jobs running after closing SSH sessions

- [Useful guide to using Tmux](https://robots.thoughtbot.com/a-tmux-crash-course)

Tmux allows you to exit an SSH session w/o killing any jobs that you initiated during an SSH session. That means, you could hop onto Emu (even with a janky internet connection), get the jobs running, and then exit the Emu connection. You can hop back onto Emu at any time to check on the status of your jobs - it'll be like you never left; all the terminal ouputs will be displayed as normal and all of that.

Here's an example of how I initiate my SSH sessions (I always use tmux to prevent problems that would be caused by a disconnection):

SSH into computer.
After connecting, type in the following (and then press Enter): ```tmux```
The terminal prompt won't look different, but you should notice a green bar is now displayed across the bottom of your terminal screen; this indicates you're in a tmux session.
Initiate any computing stuff you want.
After they're started, you can exit the SSH session in the following fashion:

1. Press ctrl-b on your keyboard (hold the control button and press the b key on the keyboard). This is a special tmux key combination to tell it to interpret the next command.
2. Press d on your keyboard. This "detaches" your tmux session. All your stuff that you started in tmux is still running in that session. When you do this, you should see a message indicating that you've detached and which tmux session you've detached from (indicated by a number; usually 0 if this is your first/only tmux session).
3. Exit your SSH session (type exit and press enter).
You've successfully initiated long-running computing junk over a janky internet connection and no longer have to worry about a bad connection interrupting your computing jobs!

To come back to your running jobs, do the following:

1. SSH into the same computer as before.
2. Attach your previous tmux session where your jobs are running (type the following and press Enter): ```tmux attach -t0```
This tells tmux to attach terminal session 0. (If you happen to have additional tmux sessions running, you can see which one you want by listing your tmux sessions with the following command: ```tmux ls```)
