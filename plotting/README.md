If you want to run these plotting scripts on your local machine using jupyter lab with the proper packages and dependencies, you can use the environment.yaml file to set up a virtual environment and jupyter kernel with the following steps:
1. Create the virtual environment with all the required packages:
   
   ```conda env create -f environment.yaml```
   
   Note: You can customize the name of this environment by editting the last line of the environment.yaml file.
3. Activate the virtual environment:
   
   ```conda activate dcmip2025```
   
5. Create a kernel using this virtual environment:
   
   ```python -m ipykernel install --user --name=dcmip2025 --display-name dcmip2025```
   
Now you should be all set up. To test this:

1. Open up JupyterLab with<br>

   ```jupyter lab```
   
3. Open/create a jupyter .ipynb file
   
4. In the top right, in between the debug and an open bug symbol is your current kernel. Click on the kernel to see other accessible kernels. The dcmip2025 kernel should now appear here. If not, try restarting the kernel and checking the list again.
