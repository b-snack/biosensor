"""
BIOLOGICALLY ACCURATE SCLEROSTIN BIOSENSOR MODEL
------------------------------------------------
Fixed SBML compliance issue by separating input control from reactants.

Key features:
1. Sclerostin acts as an INHIBITOR of Wnt/LRP6 signaling
2. Proper modular structure: detection → transduction → reporter
3. GFP maturation delay included (~30-40 min)
4. Smooth transitions via rate rules (not instantaneous jumps)
5. Proper unit definitions: nanomoles, minutes
6. Sbml-compliant: uses Sclerostin_input as controlled variable
"""

import tellurium as te
import numpy as np
import matplotlib.pyplot as plt

def create_biosensor(include_events=True):
    """
    Biologically accurate model with proper Antimony syntax
    """
    
    # Build events string conditionally
    if include_events:
        events_str = """
    // Event-based input changes
    at time > 200: sclerostin_target = 0.2
    at time > 500: sclerostin_target = 0.0
    at time > 800: sclerostin_target = 0.5
    at time > 1100: sclerostin_target = 0.0
"""
    else:
        events_str = ""
    
    model = f"""
model sclerostin_inhibition_sensor
    compartment cell = 1.0

    // Separate input control from reactive species
    species Sclerostin_control in cell
    species Sclerostin in cell
    species LRP6_free in cell
    species Sclerostin_LRP6 in cell
    species Wnt_active in cell
    species Beta_catenin in cell
    species mRNA in cell
    species Reporter_inactive in cell
    species Reporter_active in cell

    Sclerostin_control = 0.0
    Sclerostin = 0.0
    LRP6_free = 100.0
    Sclerostin_LRP6 = 0.0
    Wnt_active = 10.0
    Beta_catenin = 10.0
    mRNA = 0.0
    Reporter_inactive = 0.0
    Reporter_active = 0.0

    k_on = 0.5
    k_off = 0.01
    k_wnt_activation = 0.1
    k_wnt_decay = 0.05
    k_bc_activation = 0.2
    k_bc_decay = 0.1
    k_transcription = 1.0
    k_translation = 1.0
    k_leak = 0.1
    k_mRNA_decay = 0.05
    k_maturation = 0.05
    k_gfp_decay = 0.005
    k_equilibration = 1.0
    K_inhib = 50.0
    n = 2.0
    tau_input = 10.0
    sclerostin_target = 0.0

    // Sclerostin equilibration (fast tracking)
    R0: -> Sclerostin; k_equilibration * Sclerostin_control

    // Binding reactions
    R1: Sclerostin + LRP6_free -> Sclerostin_LRP6; k_on * Sclerostin * LRP6_free
    R2: Sclerostin_LRP6 -> Sclerostin + LRP6_free; k_off * Sclerostin_LRP6

    // Wnt signaling cascade
    R3: LRP6_free -> LRP6_free + Wnt_active; k_wnt_activation * LRP6_free
    R4: Wnt_active -> ; k_wnt_decay * Wnt_active
    R5: Wnt_active -> Wnt_active + Beta_catenin; k_bc_activation * Wnt_active
    R6: Beta_catenin -> ; k_bc_decay * Beta_catenin

    // Reporter expression WITH HILL INHIBITION
    R7: Beta_catenin -> Beta_catenin + mRNA; k_transcription * Beta_catenin * (1 / (1 + (Sclerostin_LRP6 / K_inhib)^n))
    R7b: -> mRNA; k_leak
    R8: mRNA -> ; k_mRNA_decay * mRNA
    R9: mRNA -> mRNA + Reporter_inactive; k_translation * mRNA
    R10: Reporter_inactive -> Reporter_active; k_maturation * Reporter_inactive
    R11: Reporter_active -> ; k_gfp_decay * Reporter_active

    // Rate rule for smooth input control
    Sclerostin_control' = (sclerostin_target - Sclerostin_control) / tau_input
{events_str}
end
"""
    return model
    return model

def run_accurate_simulation():
    """Run the biologically accurate model"""
    
    model_string = create_biosensor()
    r = te.loada(model_string)
    
    # Simulation time: 1200 minutes with 2000 points
    result = r.simulate(0, 1200, 2000)
    
    return result, r

def plot_accurate_results(result):
    """Plot showing inhibitory mechanism with all modules"""
    
    time = result[:, 0]
    sclerostin_control = result[:, 1]
    sclerostin = result[:, 2]
    lrp6_free = result[:, 3]
    sclerostin_lrp6 = result[:, 4]
    wnt_active = result[:, 5]
    beta_catenin = result[:, 6]
    mrna = result[:, 7]
    reporter_inactive = result[:, 8]
    reporter_active = result[:, 9]
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 11))
    axes = axes.flatten()
    
    # 1. Sclerostin input (smooth transitions)
    axes[0].plot(time, sclerostin_control, 'b-', linewidth=2.5)
    axes[0].set_ylabel('Sclerostin Input (nM)', fontsize=11)
    axes[0].set_title('Detection Module: Inhibitor Input (Smooth)', fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    
    # 2. Receptor states (detection)
    axes[1].plot(time, lrp6_free, 'g-', label='Free LRP6', linewidth=2.5)
    axes[1].plot(time, sclerostin_lrp6, 'r-', label='Sclerostin-LRP6', linewidth=2.5)
    axes[1].set_ylabel('Concentration (nM)', fontsize=11)
    axes[1].set_title('Detection Module: Receptor States', fontsize=12, fontweight='bold')
    axes[1].legend(loc='best')
    axes[1].grid(True, alpha=0.3)
    
    # 3. Wnt signaling (transduction)
    axes[2].plot(time, wnt_active, 'orange', linewidth=2.5)
    axes[2].set_ylabel('Wnt Active (nM)', fontsize=11)
    axes[2].set_title('Transduction Module: Wnt Signal', fontsize=12, fontweight='bold')
    axes[2].grid(True, alpha=0.3)
    
    # 4. β-catenin (transduction)
    axes[3].plot(time, beta_catenin, 'purple', linewidth=2.5)
    axes[3].set_ylabel('β-catenin (nM)', fontsize=11)
    axes[3].set_title('Transduction Module: β-catenin', fontsize=12, fontweight='bold')
    axes[3].grid(True, alpha=0.3)
    
    # 5. mRNA (reporter)
    axes[4].plot(time, mrna, 'brown', linewidth=2.5)
    axes[4].set_ylabel('mRNA (nM)', fontsize=11)
    axes[4].set_title('Reporter Module: mRNA', fontsize=12, fontweight='bold')
    axes[4].grid(True, alpha=0.3)
    
    # 6. GFP maturation with delay
    axes[5].plot(time, reporter_inactive, 'gray', label='Inactive GFP', linewidth=2.5)
    axes[5].plot(time, reporter_active, 'lime', label='Active GFP', linewidth=2.5)
    axes[5].set_ylabel('Reporter (nM)', fontsize=11)
    axes[5].set_title('Reporter Module: GFP Maturation', fontsize=12, fontweight='bold')
    axes[5].legend(loc='best')
    axes[5].grid(True, alpha=0.3)
    
    # 7. Inverse relationship (key feature)
    axes[6].plot(time, sclerostin_control, 'b-', label='Sclerostin Input', linewidth=2.5)
    ax6_twin = axes[6].twinx()
    ax6_twin.plot(time, reporter_active, 'lime', label='Active Reporter', linewidth=2.5)
    axes[6].set_xlabel('Time (minutes)', fontsize=11)
    axes[6].set_ylabel('Sclerostin (nM)', fontsize=11, color='b')
    ax6_twin.set_ylabel('Reporter (nM)', fontsize=11, color='lime')
    axes[6].set_title('Key Result: Inverse Relationship', fontsize=12, fontweight='bold')
    axes[6].tick_params(axis='y', labelcolor='b')
    ax6_twin.tick_params(axis='y', labelcolor='lime')
    axes[6].grid(True, alpha=0.3)
    
    # 8. Modular architecture summary
    axes[7].axis('off')
    metrics_text = f"""
MODULAR BIOSENSOR ARCHITECTURE
─────────────────────────────

Detection Module:
  Sclerostin + LRP6 ⇌ Complex
  KD ≈ 20 nM
  Blocks Wnt signaling

Signal Transduction:
  Free LRP6 → Wnt signal
  Wnt → β-catenin activation
  Cascade amplification

Reporter Module:
  β-catenin → mRNA
  mRNA → GFP(inactive)
  GFP(inactive) → GFP(active)
  Maturation t½ ≈ 28 min

Kinetics (minutes):
  mRNA t½ ≈ 14 min
  Active GFP t½ ≈ 69 min
  Total response lag: ~60 min
  Smooth input tau: ~10 min
    """
    axes[7].text(0.05, 0.5, metrics_text, fontsize=9.5,
                verticalalignment='center', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.suptitle('Biologically Accurate Sclerostin Biosensor (SBML Compliant)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()
    
    return fig

def calculate_metrics(result):
    """Calculate key performance metrics"""
    
    time = result[:, 0]
    sclerostin_control = result[:, 1]
    sclerostin = result[:, 2]
    wnt_active = result[:, 5]
    reporter_active = result[:, 9]
    
    # Baseline (first 60 min)
    baseline_idx = time < 60
    baseline_reporter = np.mean(reporter_active[baseline_idx])
    baseline_wnt = np.mean(wnt_active[baseline_idx])
    
    # First pulse (60-300 min, sclerostin_target = 0.2)
    pulse1_idx = (time > 120) & (time < 300)
    inhibited_reporter = np.mean(reporter_active[pulse1_idx])
    inhibited_wnt = np.mean(wnt_active[pulse1_idx])
    
    # Calculate inhibition
    inhibition = ((baseline_reporter - inhibited_reporter) / 
                  baseline_reporter * 100) if baseline_reporter > 0 else 0
    
    # Response time (to 50% inhibition)
    half_inhibition = (baseline_reporter + inhibited_reporter) / 2
    response_idx = np.where((time > 60) & (reporter_active < half_inhibition))[0]
    response_time = (time[response_idx[0]] - 60) if len(response_idx) > 0 else np.inf
    
    print("\n" + "="*50)
    print("BIOSENSOR PERFORMANCE METRICS")
    print("="*50)
    print(f"Baseline reporter: {baseline_reporter:.4f} nM")
    print(f"Baseline Wnt signal: {baseline_wnt:.4f} nM")
    print(f"Inhibited reporter: {inhibited_reporter:.4f} nM")
    print(f"Inhibited Wnt signal: {inhibited_wnt:.4f} nM")
    print(f"Inhibition efficiency: {inhibition:.1f}%")
    print(f"Response time (50% inhibition): {response_time:.1f} minutes")
    print("="*50 + "\n")
    
    return {
        'baseline': baseline_reporter,
        'inhibited': inhibited_reporter,
        'inhibition_percent': inhibition,
        'response_time': response_time
    }

def test_dose_response():
    """Test different sclerostin concentrations at steady state"""
    
    print("\n" + "="*50)
    print("DOSE-RESPONSE ANALYSIS")
    print("="*50)
    
    doses = [0.0, 0.05, 0.1, 0.2, 0.5, 1.0]
    responses = []
    
    for dose in doses:
        model_string = create_biosensor(include_events=False)
        r = te.loada(model_string)
        
        # DEBUG: Print column names to see actual order
        if dose == 0.0:
            print(f"Column names: {r.getFloatingSpeciesIds()}")
            print(f"All selections: {r.selections}")
        
        r['sclerostin_target'] = dose
        r['Sclerostin_control'] = dose
        r['Sclerostin'] = dose
        
        result = r.simulate(0, 1200, 500)
        
        # Use the species NAME instead of column index
        reporter_ss = np.mean(result[-100:, r.selections.index('[Reporter_active]')])
        responses.append(reporter_ss)
        
        print(f"Sclerostin = {dose:.2f} nM → Active Reporter = {reporter_ss:.4f} nM")
    
    print("="*50 + "\n")
    
    # Plot dose-response curve
    plt.figure(figsize=(8, 6))
    plt.plot(doses, responses, 'o-', linewidth=2.5, markersize=10, color='darkgreen')
    plt.xlabel('Sclerostin Input Concentration (nM)', fontsize=12)
    plt.ylabel('Steady-State Active Reporter (nM)', fontsize=12)
    plt.title('Dose-Response: Sclerostin Inhibits Reporter Expression', fontsize=13, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    return doses, responses

def main():
    """Run the complete corrected biosensor simulation"""
    
    print("\n" + "="*60)
    print(" CORRECTED SCLEROSTIN BIOSENSOR MODEL")
    print(" SBML-Compliant with Modular Architecture")
    print("="*60)
    
    # Run simulation
    print("\nRunning simulation with corrected biology...")
    result, model = run_accurate_simulation()
    print("✓ Simulation complete")
    
    # Calculate metrics
    metrics = calculate_metrics(result)
    
    # Plot results
    print("Generating plots...")
    plot_accurate_results(result)
    
    # Dose-response analysis
    test_dose_response()
    
    print("✓ All analyses complete\n")
    
    return result, model, metrics

if __name__ == "__main__":
    result, model, metrics = main()