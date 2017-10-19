## Tips and tricks for connecting to Roberts Lab computers remotely.

### Copy Files from/to Locations Using SSH (Terminal)
 (NOTE: This will not work natively in Windows, as Windows doesn't have an SSH client. Download and install Cygwin or OpenSSH to gain SSH functionality in Windows).


 Copy something from this system to some other system:
 `scp /path/to/local/file username@hostname:/path/to/remote/file`


 Copy something from some system to some other system:
 `scp username1@hostname1:/path/to/file username2@hostname2:/path/to/other/file`


 Copy something from another system to this system:
 `scp username@hostname:/path/to/remote/file /path/to/local/file`

 ---


### Set Up SSH Keys for More Secure SSH to Hummingbird(Terminal)
 For any SSH connection to Humminbird (or any server), SSH keys should be used.

 _Instructions are for Macintosh OS X_
 1. Generate your SSH keys.
 In Terminal, paste the following:
 `ssh-keygen -t rsa`    
 Press "Enter" at both prompts for a passphrase.
 Note: You can certainly enter a passphrase if desired, but the passphrase is only needed in instances where someone has obtained control of your laptop/desktop that you use to connect to Hummingbird. By not entering a passphrase, you obtain the luxury of being able to use SSH without the need to type a password.

 2. Obtain a script to simplify copying your SSH key to the server.
 In Terminal, paste the following:
 ```
 sudo curl https://raw.github.com/beautifulcode/ssh-copy-id-for-OSX/master/ssh-copy-id.sh -o /usr/local/bin/ssh-copy-id -k
 ```

 3. Change permissions on the "ssh-copy-id" script to make it executable.
 In Terminal, paste the following:
 ```
 sudo chmod +x /usr/local/bin/ssh-copy-id
 ```

 4. Copy your SSH key to the server using the "ssh-copy-id" script.
 In Terminal, paste the following:
 ```
 ssh-copy-id username@Hummingbird's.IP.address
 ```
 You'll be prompted to enter your password for your Hummingbird user account.

 5. Test it out.
 In Terminal, connect to Humminbird using SSH

---

### Tmux to keep jobs running after closing SSH sessions

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
