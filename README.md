# 20251222_Zesheng_Xuanwu_brain
MFA for brain metabolism


## Getting started

This script could also be executed as a raw Python project. Make sure Python 3.8 and all required packages are correctly
installed. First switch to a target directory and download the source code:

```shell script
git clone https://github.com/cmplab-cimr/2025_neurodegenerative_metabolism
```

Switch to the source direct, add PYTHONPATH environment and run the `main.py`:

```shell script
cd 2025_neurodegenerative_metabolism
export PYTHONPATH=$PYTHONPATH:`pwd`
python main.py
```

You could try multiple different arguments according to help information. Different mode is support for multiple functions, including:


```shell script
# Experimental data plotting:
python main.py computation experimental_mfa raw_experimental_data_plotting

# Flux analysis:
python main.py computation experimental_mfa flux_analysis -p 20

# Process results of flux analysis:
python main.py computation experimental_mfa result_process
```

This code support running in Docker. The script `docker_commands.sh` is provided for convenient operations, including:

```shell script
# Start container
bash docker_commands.sh -p 20 start

# View container logs
bash docker_commands.sh log

# Stop container
bash docker_commands.sh stop

# Process results of flux analysis
bash docker_commands.sh result

# Build image
bash docker_commands.sh build
```
