

You need to use `tmux` to run a long-running process, otherwise if your ssh connection drops, the program will not continue.

after you log in to the computer, type
`tmux`

You should see a bar at the bottom of the screen, 

Now, you can run the program as you normally would. The difference will be, if by any chance you drop off ssh, your program will keep on running, and you can see it after reconnecting.

Say if you've closed the window that you are running the program in, or got disconnected, connect again and type

`tmux a`

This will reconnect to the previous session, and you'll be able to see the program running.

If you want to see if any previous session is running, you can check by typing:

`tmux list-sessions`

You can close the program as you would by pressing `Ctrl+C`