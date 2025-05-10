# Setting up VSCode Remote

Normally, VSCode only allows you to edit/view files locally. However, you can configure VSCode-Remote to work with files stored in virtual machines/servers.

1. Install the **Remote-SSH** extension on VSCode.

2. Ensure you have a SSH key, which you can access by checking `cat ~/.ssh/id_ed25519` on the server. If the file doesn't exist, you can create it by:

```bash
ssh-keygen -t ed25519
```

3. Open your Command Palette in VSCode by typing `CTRL + SHIFT + P`. This should open a prompt. Enter `Connect to Host` which should have `Remote-SSH` next to it (if it doesn't, verify that you have **Remote-SSH** installed).

4. Enter your connection information in the form `user@hostname`. For example, for Cedar, it would be:

```
[username]@cedar.alliancecan.ca
```

5. Open the **Remote Explorer** (Computer with a little connection icon in the corner). You will see a list of options under **Remote**. Click on the your hostname listed under **SSH**.

6. This opens a new instance of VSCode. Enter your password as prompted to in the search bar on top.

7. Navigate to the **Terminal** window in the bottom of VSCode to complete the two-factor authentication. This is kind of finnicky, so if it is empty, click on the "details" link in the prompt that opens in the bottom right. It will say "Setting up SSH Host [hostname]: (details) Initializing VSCode Server." in this prompt. 

You should see the option to trigger 2FA, just as you would when doing `ssh` into the server from your terminal. You may be prompted to enter your password in VSCode a second time.

8. If all is well, you should now be able to select **Open Folder** in the left sidebar. You may need to recomplete the user login/authentication step after selecting the folder you want to open.

9. Note that opening files can sometimes not work. You can use `CTRL + p` and type the name of the file to open it if clicking on the file in the left sidebar doesn't work.

10. To access the server shell from within VSCode, you may need to open a new shell (**Terminal** window in VSCode). Click the little "+" icon in the top right of this window to add a new terminal.
