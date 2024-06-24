import os
import shutil
import subprocess
import re

num_iterations = 16  # Choose number of iteration wanted
initial_temp = 0.45  # Initial temperature as set in input.dat

for i in range(num_iterations):

    try:
        subprocess.run(["./simulator.exe"], check=True)
        print("Iteration number " + str(i+1) + " successfully executed!")
    except subprocess.CalledProcessError as e:
        print("Error during simulator.exe execution:", e)
        exit(1)

    
    # Modify input.dat after the first execution
    with open("../INPUT/input.dat", "r") as input_file:
        lines = input_file.readlines()
    if i == 0:
        with open("../INPUT/input.dat", "w") as input_file:
            for line in lines:
                if "RESTART" in line:
                    line = "RESTART                1\n"  # RESTART becomes 1 after first iteration
                input_file.write(line)

    # Move and rename .dat files in ../OUTPUT
    output_folder = "../OUTPUT"
    result_folder = "../RESULTS/Esercitazione_06/Iteration_data"
    for file in os.listdir(output_folder):
        if file.endswith(".dat"):
            file_name = os.path.splitext(file)[0]
            shutil.move(os.path.join(output_folder, file), os.path.join(result_folder, f"{file_name}_{i+1}.dat"))
            file_numbered_path = os.path.join(result_folder, f"{file_name}_{i+1}.dat")
            file_unnumbered_path = os.path.join(result_folder, f"{file_name}.dat")
            # Create a file for all last lines only after the firts iteration
            if i == 0 and not os.path.exists(file_unnumbered_path):
                with open(file_unnumbered_path, 'w') as f:
                    f.write("#      TEMP:         ACTUAL:             AVE:             ERROR:\n")
            with open(file_numbered_path, 'r') as filein:
                lines = filein.readlines()
                last_line = lines[-1]
                first_item = last_line.split()[0]
                temp = "{:.2f}".format(initial_temp + i * 0.1)
                new_last_line = last_line.replace(first_item, str(temp), 1)
            with open(file_unnumbered_path, 'a') as f:
                f.write(new_last_line)

    # Modify TEMP after every iteration
    new_temp = "{:.2f}".format(initial_temp + (i+1) * 0.1) 
    with open("../INPUT/input.dat", "r") as input_file:
        lines = input_file.readlines()
    with open("../INPUT/input.dat", "w") as input_file:
        for line in lines:
            if "TEMP" in line:
                line = f"TEMP                   {new_temp}\n"
            input_file.write(line)

    # Move config.spin in ../INPUT/CONFIG
    input_config_folder = "../INPUT/CONFIG"
    output_config_folder = "../OUTPUT/CONFIG"
    for file in os.listdir(output_config_folder):
        if file == "config.spin":
            shutil.move(os.path.join(output_config_folder, file), os.path.join(input_config_folder, "config.spin"))

# Reset RESTART to 0 and TEMP to 0.45 after all iterations
with open("../INPUT/input.dat", "r") as input_file:
    lines = input_file.readlines()
modified_lines = []
for line in lines:
    if "RESTART" in line:
        line = "RESTART                0\n"  # Reset RESTART to 0
    if "TEMP" in line:
        line = "TEMP                   0.45\n"  # Reset TEMP to 0.45
    modified_lines.append(line)
with open("../INPUT/input.dat", "w") as input_file:
    input_file.writelines(modified_lines)


print("Iterations completed.")
