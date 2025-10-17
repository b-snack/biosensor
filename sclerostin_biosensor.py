"""
BASIC SCLEROSTIN BIOSENSOR MODEL
---------------------------------
Minimal version with just the essential components:
1. Sclerostin detection
2. Signal generation
3. Reporter output
"""

import tellurium as te
import numpy as np
import matplotlib.pyplot as plt

# Create the basic model
def create_basic_model():
    """
    Minimal biosensor with:
    - Sclerostin input
    - Receptor binding
    - Reporter output
    """
    
    model = """
    model basic_sclerostin_sensor
        
        # Compartment
        compartment cell = 1.0;
        
        # Species
        species Sclerostin in cell;    # Input molecule we want to detect
        species Receptor in cell;       # Sensor that binds sclerostin
        species Complex in cell;        # Bound receptor-sclerostin
        species Reporter in cell;       # Output signal (e.g., GFP)
        
        # Initial values
        Sclerostin = 0.0;   # Start with no sclerostin
        Receptor = 10.0;    # Available receptors
        Complex = 0.0;      # No bound complex initially
        Reporter = 0.0;     # No output initially
        
        # Parameters
        k_bind = 1.0;       # Binding rate
        k_unbind = 0.1;     # Unbinding rate  
        k_produce = 2.0;    # Reporter production rate
        k_degrade = 0.5;    # Reporter degradation rate
        
        # Reactions
        R1: Sclerostin + Receptor -> Complex;  k_bind * Sclerostin * Receptor;
        R2: Complex -> Sclerostin + Receptor;  k_unbind * Complex;
        R3: Complex -> Complex + Reporter;     k_produce * Complex;
        R4: Reporter -> ;                       k_degrade * Reporter;
        
        # Test events - sclerostin spikes
        at (time >= 10): Sclerostin = 2.0;   # Add sclerostin
        at (time >= 30): Sclerostin = 0.0;   # Remove sclerostin
        at (time >= 50): Sclerostin = 4.0;   # Higher dose
        at (time >= 70): Sclerostin = 0.0;   # Remove again
        
    end
    """
    return model

# Run simulation
def run_basic_simulation():
    """Run the basic model and return results"""
    
    # Load model
    model_string = create_basic_model()
    r = te.loada(model_string)
    
    # Run simulation
    result = r.simulate(0, 100, 500)
    
    return result, r

# Plot results
def plot_basic_results(result):
    """Simple plot of the key species"""
    
    time = result[:, 0]
    sclerostin = result[:, 1]
    receptor = result[:, 2]
    complex = result[:, 3]
    reporter = result[:, 4]
    
    fig, axes = plt.subplots(4, 1, figsize=(10, 10))
    
    # Sclerostin input
    axes[0].plot(time, sclerostin, 'b-', linewidth=2)
    axes[0].set_ylabel('Sclerostin')
    axes[0].set_title('Input Signal')
    axes[0].grid(True, alpha=0.3)
    
    # Free receptor
    axes[1].plot(time, receptor, 'g-', linewidth=2)
    axes[1].set_ylabel('Free Receptor')
    axes[1].set_title('Available Receptors')
    axes[1].grid(True, alpha=0.3)
    
    # Bound complex
    axes[2].plot(time, complex, 'orange', linewidth=2)
    axes[2].set_ylabel('Bound Complex')
    axes[2].set_title('Receptor-Sclerostin Complex')
    axes[2].grid(True, alpha=0.3)
    
    # Reporter output
    axes[3].plot(time, reporter, 'r-', linewidth=2)
    axes[3].set_ylabel('Reporter')
    axes[3].set_xlabel('Time')
    axes[3].set_title('Output Signal')
    axes[3].grid(True, alpha=0.3)
    
    plt.suptitle('Basic Biosensor Response', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()
    
    return fig

# Calculate simple metrics
def calculate_basic_metrics(result):
    """Calculate basic performance metrics"""
    
    time = result[:, 0]
    sclerostin = result[:, 1]
    reporter = result[:, 4]
    
    # Find maximum reporter signal
    max_reporter = np.max(reporter)
    
    # Find baseline (before first spike)
    baseline = np.mean(reporter[time < 10])
    
    # Find signal during first spike
    spike1_signal = np.max(reporter[(time >= 10) & (time <= 30)])
    
    # Calculate fold change
    if baseline > 0:
        fold_change = spike1_signal / baseline
    else:
        fold_change = spike1_signal
    
    print("\n=== BASIC METRICS ===")
    print(f"Baseline signal: {baseline:.3f}")
    print(f"Max reporter signal: {max_reporter:.3f}")
    print(f"First spike signal: {spike1_signal:.3f}")
    print(f"Fold change: {fold_change:.1f}x")
    print("=" * 20)

# Main execution
def main():
    """Run the basic biosensor simulation"""
    
    print("\n" + "="*50)
    print(" BASIC SCLEROSTIN BIOSENSOR")
    print("="*50 + "\n")
    
    # Run simulation
    print("Running simulation...")
    result, model = run_basic_simulation()
    print("✓ Simulation complete")
    
    # Show what species we have
    print("\nModel species:")
    for i, species in enumerate(model.getFloatingSpeciesIds()):
        print(f"  {i}: {species}")
    
    # Calculate metrics
    calculate_basic_metrics(result)
    
    # Plot results
    print("\nGenerating plot...")
    plot_basic_results(result)
    
    return result, model

# Parameter testing function
def test_parameter(param_name, test_values):
    """Test how changing a parameter affects the output"""
    
    print(f"\n=== Testing {param_name} ===")
    
    for value in test_values:
        # Create and load model
        model_string = create_basic_model()
        r = te.loada(model_string)
        
        # Set parameter
        r[param_name] = value
        
        # Run simulation
        result = r.simulate(0, 100, 500)
        
        # Get max reporter value
        max_reporter = np.max(result[:, 4])
        
        print(f"{param_name} = {value:.2f} → Max signal = {max_reporter:.3f}")

if __name__ == "__main__":
    # Run main simulation
    result, model = main()
    
    # Optional: Test parameters
    print("\n" + "="*50)
    print(" PARAMETER TESTING")
    print("="*50)
    
    # Test binding affinity
    test_parameter('k_bind', [0.5, 1.0, 2.0, 5.0])
    
    # Test production rate
    test_parameter('k_produce', [1.0, 2.0, 5.0, 10.0])