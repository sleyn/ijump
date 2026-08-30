## Running the Simulation

This guide details how to create test data for ijump

**Prerequisites:**

* Linux system

**Steps:**

1. Navigate to the simulation directory:

```bash
cd ijump/simulation
```

2. Run the script to create test input data:

```bash
bash create_test_input_ijump.sh
```

# Notes

1. The alignment.sh script (line 17) is configured to use 6 threads (-t 6).
2. This script is currently only functional on Linux due to its reliance on apt-get install.