import subprocess
import time
import matplotlib.pyplot as plt

def run_benchmark(mol_sim_path, input_file, benchmark_option="", newton_option="", flag_option="", num_iterations=10):
    execution_times = []
    
    for _ in range(num_iterations):
        start_time = time.time()
        
        # Execute Molsim.cpp with the given options
        result = subprocess.run(
            [mol_sim_path, input_file, benchmark_option, newton_option, flag_option], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )
        
        # Measure execution time
        end_time = time.time()
        execution_times.append(end_time - start_time)
        
        # Check for errors
        if result.returncode != 0:
            print("Error:", result.stderr.decode())
            return None
    
    # Return the average execution time over all iterations
    avg_time = sum(execution_times) / len(execution_times)
    return avg_time

def plot_benchmark(time_with_10_l, time_without_10_l, time_with_100_l, time_without_100_l, time_with_1000_l, time_without_1000_l):
    labels = ["O 10", "X 10", "O 100", "X 100", "O 1000", "X 1000"]
    
    times = [time_with_10_l, time_without_10_l, time_with_100_l, time_without_100_l, time_with_1000_l, time_without_1000_l]

    # Set font size for titles and labels
    plt.rcParams.update({'font.size': 10})

    plt.bar(labels, times, color=['blue', 'red', 'blue', 'red', 'blue', 'red'])
    plt.title('Average Execution Time for Molsim.cpp with and without Newton\'s 3rd law(LennardJonesForce)', fontsize=8)
    plt.ylabel('Average Execution Time (seconds)', fontsize=10)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(fontsize=9)
    plt.grid(True)
    plt.show()

def main():
    # Path
    mol_sim_path = "./build/MolSim"
    input_file = "./input/eingabe-sonne.txt"

    # Run benchmark with Newton's 3rd law + -l flag for different iterations
    avg_time_with_10_l = run_benchmark(mol_sim_path, input_file, benchmark_option="-b", newton_option="", flag_option="-l", num_iterations=10)
    if avg_time_with_10_l is None:
        return
    avg_time_without_10_l = run_benchmark(mol_sim_path, input_file, benchmark_option="-b", newton_option="-n", flag_option="-l", num_iterations=10)
    if avg_time_without_10_l is None:
        return

    avg_time_with_100_l = run_benchmark(mol_sim_path, input_file, benchmark_option="-b", newton_option="", flag_option="-l", num_iterations=100)
    if avg_time_with_100_l is None:
        return
    avg_time_without_100_l = run_benchmark(mol_sim_path, input_file, benchmark_option="-b", newton_option="-n", flag_option="-l", num_iterations=100)
    if avg_time_without_100_l is None:
        return

    avg_time_with_1000_l = run_benchmark(mol_sim_path, input_file, benchmark_option="-b", newton_option="", flag_option="-l", num_iterations=1000)
    if avg_time_with_1000_l is None:
        return
    avg_time_without_1000_l = run_benchmark(mol_sim_path, input_file, benchmark_option="-b", newton_option="-n", flag_option="-l", num_iterations=1000)
    if avg_time_without_1000_l is None:
        return

    # Plot average execution times for comparison across 10, 100, and 1000 iterations with -l flag
    plot_benchmark(
        avg_time_with_10_l, avg_time_without_10_l, avg_time_with_100_l, avg_time_without_100_l, avg_time_with_1000_l, avg_time_without_1000_l
    )

if __name__ == "__main__":
    main()
