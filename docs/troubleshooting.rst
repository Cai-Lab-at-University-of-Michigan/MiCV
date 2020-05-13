Troubleshooting
================================

This file describes some of the common issues and problem solving
techniques that are relevant to scRNA-seq analysis in MiCV. It will grow over time as we uncover and try to address new issues!

----------

The interface appears to be frozen
**********************************

- Did you just upload data? If so, intial upload and data parsing might take a while - check your title bar to see if it says "updating", or check a network monitor to see if an upload is still in progress.
- Check the status bar on the right-hand side of the screen, as well as the title bar for the tab of your web-browser. Is the status bar updating, or does the title bar say "updating"? If so, there is likely a long-running calculation going on the background - wait a few more minutes to see if things update.
- Datasets with lots of cells make everything, including plot updates, take longer - we're working to speed this up, but in the mean time, please be patient or consider sub-sampling your cells. 
- The interface may have frozen due to a bug in the code (often happens if an algorithm fails to converge, like in pseudotime calculations) - you'll have to refresh your tab and start over from the beginning for now! Please let us know what happened just before the crash so we can look into what went wrong and how to prevent that error in the future.

Data upload fails completely
****************************

- We have found that when uploading ``h5ad`` files to MiCV, Google Chrome tends to throw a security-related error. The root cause of this is not known to us yet, but we know that Mozilla Firefox seems not to have this issue. As such, we recommend users interact with MiCV through Firefox for the time being.