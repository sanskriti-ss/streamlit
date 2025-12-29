# Guideline for downloading circos files:
Go to the Circos tabs on https://bacdive.streamlit.app/ and make your selections.

Note that the data available on the streamlit is from the latest data download done (data available in readme). If you want to use a more recent BacDive database date, you must run it yourself in https://github.com/AshishMahabal/antibiotic

# Instructions for running on circos:
1) Make sure you have Circos properly installed onto your device, do the test run as the documentation suggests, etc.
2) Enter the specifications of what you want to visualize on the Streamlit by using the dropdown and sliders.
3) Download the generated files.
4) Navigate to Circos using the terminal, and then the ‘bin’. Copy and paste the text from the downloaded file with the suffix “chromosome.txt" into the karyotype file. Save.
5) Make sure that you’re using the right file name (for the downloaded file with the ‘links’). Move it into this folder. 
    (base) PS C:\Users\usern\circos-0.69-9\bin> ./circos -conf /path/to/top_fifty_circos2.conf -noparanoid
Note that you can put specific things in the conf file, but you don’t want them to do weird things with the settings we are downloading. I recommend you just use mine (will link here soon). Every time you want to run something, you’ll need to change the name of the link file. So if you downloaded top_genera_res_1strain_top_50_Triples_.txt from Streamlit, change the file name in this line in the config file: file = top_fifty_link_data.txt to the proper file name, and make sure your downloaded file is in this folder!!!!

6) If you downloaded secondary files, you need to copy and paste the links into the bottom of the primary file (don’t overwrite any of the links in there) as you can only give circos one links file; so you must manually combine them.

7) Run it. It might take some time the first time, depends on how many links you have and all.

