import subprocess
import time
import matplotlib.pyplot as plt

def run_benchmark(mol_sim_path, input_file, benchmark_option="", newton_option=""):
    execution_times = []
    
    for _ in range(10):
        start_time = time.time()
        
        # Execute Molsim.cpp with the given options
        result = subprocess.run(
            [mol_sim_path, input_file, benchmark_option, newton_option], 
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

def plot_benchmark(time_with, time_without):
    # Plot graph comparing both modes
    labels = ["With Newton's 3rd law", "Without Newton's 3rd law"]
    times = [time_with, time_without]
    
    plt.bar(labels, times, color=['blue', 'red'])
    plt.title('Benchmark Comparison for MolSim')
    plt.ylabel('Average Execution Time (seconds)')
    plt.grid(True)
    plt.show()

def main():
    # Path
    mol_sim_path = "./build/MolSim"
    input_file = "./input/eingabe-sonne.txt"

    # Run benchmark with Newton's 3rd law (without -n option)
    avg_time_with = run_benchmark(mol_sim_path, input_file, benchmark_option="-b", newton_option="")
    if avg_time_with is None:
        return
    
    # Run benchmark without Newton's 3rd law (with -n option)
    avg_time_without = run_benchmark(mol_sim_path, input_file, benchmark_option="-b", newton_option="-n")
    if avg_time_without is None:
        return
    
    # Plot average execution times for comparison
    plot_benchmark(avg_time_with, avg_time_without)

if __name__ == "__main__":
    main()
