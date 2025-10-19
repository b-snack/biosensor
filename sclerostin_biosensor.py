"""
SCLEROSTIN BIOSENSOR ENVIRONMENTAL SIMULATION (STABLE VERSION)
---------------------------------------------------------------
Fixed numerical stability issues:
- Rescaled environmental species to nM range
- Smoothed temperature/pH transitions
- Reduced stiffness in resource dynamics
- Added absolute/relative tolerances
"""


import tellurium as te
import numpy as np
import matplotlib.pyplot as plt


def create_biosensor_with_environment():
    """
    Numerically stable biosensor with environmental factors
    """
    
    model = """
model biosensor_in_environment8
    compartment cell = 1.0
    
    // ============ BIOSENSOR SPECIES ============
    species Sclerostin in cell
    species LRP6_free in cell
    species Sclerostin_LRP6 in cell
    species Wnt_active in cell
    species Beta_catenin in cell
    species mRNA in cell
    species Reporter_inactive in cell
    species Reporter_active in cell
    
    // ============ ENVIRONMENTAL SPECIES (RESCALED) ============
    species Glucose in cell          // Rescaled to nM
    species ATP in cell              // Rescaled to nM
    species AminoAcids in cell       
    species ROS in cell              
    species StressProtein in cell    // Heat shock proteins
    
    // ============ ENVIRONMENTAL CONTROLS (smooth transitions) ============
    species Temperature in cell      // Dynamic temperature
    species pH_level in cell         // Dynamic pH
    
    // ============ INITIAL CONDITIONS ============
    // Biosensor
    Sclerostin = 0.0
    LRP6_free = 100.0
    Sclerostin_LRP6 = 0.0
    Wnt_active = 0.0
    Beta_catenin = 0.0
    mRNA = 0.0
    Reporter_inactive = 0.0
    Reporter_active = 0.0
    
    // Environment (RESCALED for numerical stability)
    Glucose = 100.0          // nM (same scale as biosensor)
    ATP = 5000.0              // nM (rescaled from 1000)
    AminoAcids = 50.0        // nM
    ROS = 0.01                // nM
    StressProtein = 1.0      // nM
    
    Temperature = 37.0       // °C
    pH_level = 7.4           // pH units
    
    // ============ ENVIRONMENTAL PARAMETERS ============
    // Target values (for smooth transitions)
    Temperature_target = 37.0
    pH_target = 7.4
    tau_env = 20.0           // Time constant for smooth changes (minutes)
    
    Oxygen = 21.0            // % O2 (constant for now)
    
    // Temperature correction (Q10 model)
    T_ref = 37.0
    Q10 = 2.0
    temp_factor = Q10^((Temperature - T_ref)/10.0)
    
    // pH correction (Gaussian around optimal)
    pH_optimal = 7.4
    pH_sigma = 0.8
    ph_factor = exp(-((pH_level - pH_optimal)^2)/(2*pH_sigma^2))
    
    // Oxygen factor
    O2_Km = 5.0
    oxygen_factor = Oxygen / (O2_Km + Oxygen)
    
    // ============ BIOSENSOR PARAMETERS (corrected) ============
    k_on_base = 0.5
    k_off_base = 0.01
    k_wnt_activation_base = 0.05
    k_wnt_decay = 0.05
    k_bc_activation = 0.2
    k_bc_decay = 0.1
    k_transcription_base = 0.08
    k_translation_base = 0.06
    k_leak = 0.0001
    k_mRNA_decay = 0.1
    k_maturation = 0.05
    k_gfp_decay_base = 0.005
    ATP_transcription_cost_factor = 0.02
    ATP_translation_cost_factor = 0.05
    AA_translation_cost_factor = 0.02
    
    // Apply environmental corrections
    k_on = k_on_base * temp_factor * ph_factor
    k_off = k_off_base * temp_factor
    k_wnt_activation = k_wnt_activation_base * temp_factor * ph_factor
    k_transcription = k_transcription_base * temp_factor * ph_factor
    k_translation = k_translation_base * temp_factor * stress_inhibition
    k_gfp_decay = k_gfp_decay_base * (1 + 0.5 * ROS)  // ROS accelerates decay
    
    // ============ ENVIRONMENTAL PARAMETERS ============
    // External pools (constant)
    Glucose_external = 100.0
    AA_external = 50.0
    
    // Import/export rates
    k_glucose_import = 0.2
    k_aa_import = 0.1
    
    // Metabolism (REDUCED for stability)
    k_glucose_metabolism = 0.02
    k_atp_production = 2.0
    k_atp_basal = 0.5
    k_atp_synthesis_basal = 0.1
    k_aa_consumption = 0.005
    
    // Stress dynamics
    k_ros_production = 0.0005
    k_ros_scavenge = 0.08
    k_stress_induction = 0.02
    k_stress_decay = 0.01
    
    // Resource availability factors (Michaelis-Menten)
    ATP_Km = 1000.0
    AA_Km = 10.0
    atp_availability = ATP / (ATP + ATP_Km)
    aa_availability = AminoAcids / (AminoAcids + AA_Km)
    stress_inhibition = 1.0 / (1.0 + 0.5 * StressProtein)
    
    // ============ BIOSENSOR REACTIONS ============
    // Binding
    R1: Sclerostin + LRP6_free -> Sclerostin_LRP6; k_on * Sclerostin * LRP6_free
    R2: Sclerostin_LRP6 -> Sclerostin + LRP6_free; k_off * Sclerostin_LRP6
    
    // Signaling (FREE LRP6 activates Wnt) - CATALYTIC REACTION
    R3: LRP6_free -> LRP6_free + Wnt_active; k_wnt_activation * LRP6_free
    R4: Wnt_active -> ; k_wnt_decay * Wnt_active
    R5: Wnt_active -> Wnt_active + Beta_catenin; k_bc_activation * Wnt_active
    R6: Beta_catenin -> ; k_bc_decay * Beta_catenin
    
    // Transcription (ATP-limited)
    R7: Beta_catenin -> Beta_catenin + mRNA; k_transcription * Beta_catenin * atp_availability
    R7b: -> mRNA; k_leak
    R8: mRNA -> ; k_mRNA_decay * mRNA
    
    // Translation (ATP + AA limited)
    R9: mRNA -> mRNA + Reporter_inactive; k_translation * mRNA * atp_availability * aa_availability
    
    // Maturation
    R10: Reporter_inactive -> Reporter_active; k_maturation * Reporter_inactive
    R11: Reporter_active -> ; k_gfp_decay * Reporter_active
    
    // ============ ENVIRONMENTAL REACTIONS (SIMPLIFIED) ============
    // Glucose dynamics
    E1: -> Glucose; k_glucose_import * (Glucose_external - Glucose)
    E2: Glucose -> ; k_glucose_metabolism * Glucose
    
    // ATP production (simplified, no depletion reactions)
    E3: -> ATP; k_atp_production * Glucose * oxygen_factor + k_atp_synthesis_basal
    E4: ATP -> ; k_atp_basal * ATP
    E5: ATP -> ; ATP_transcription_cost_factor * k_transcription * Beta_catenin * atp_availability
    E6: ATP -> ; ATP_translation_cost_factor * k_translation * mRNA * atp_availability
    
    // Amino acid dynamics
    E7: -> AminoAcids; k_aa_import * (AA_external - AminoAcids)
    E8: AminoAcids -> ; k_aa_consumption * AminoAcids
    E9: AminoAcids -> ; AA_translation_cost_factor * k_translation * mRNA * aa_availability
    
    // Stress response
    E10: -> ROS; k_ros_production * Glucose
    E11: ROS -> ; k_ros_scavenge * ROS * (1 + StressProtein)  // Stress proteins help
    E12: ROS -> ROS + StressProtein; k_stress_induction * ROS
    E13: StressProtein -> ; k_stress_decay * StressProtein
    
    // ============ RATE RULES (smooth environmental changes) ============
    Temperature' = (Temperature_target - Temperature) / tau_env
    pH_level' = (pH_target - pH_level) / tau_env
    
    // ============ EVENTS (change targets, not direct jumps) ============
    at time == 100: Sclerostin = 0.05
    at time == 300: Temperature_target = 30.0
    at time == 400: Temperature_target = 37.0
    at time == 500: Sclerostin = 0.5
    at time == 700: pH_target = 6.8
    at time == 800: pH_target = 7.4
    at time == 900: Glucose_external = 30.0
    at time == 1000: Glucose_external = 100.0
    at time == 1100: ROS = 1.5
    at time == 1110: ROS = 2.5
    at time == 1120: ROS = 3.0
    at time == 1200: Sclerostin = 0.0
        
end
"""
    return model


def run_environmental_simulation():
    """Run with increased tolerances for stability"""
    
    print("\n" + "="*70)
    print(" SCLEROSTIN BIOSENSOR ENVIRONMENTAL SIMULATION")
    print(" Realistic cell culture and patient sample conditions")
    print("="*70)
    
    model_string = create_biosensor_with_environment()
    r = te.loada(model_string)
    
    # Configure integrator for stability
    r.integrator = 'cvode'
    r.integrator.absolute_tolerance = 1e-10
    r.integrator.relative_tolerance = 1e-6
    r.integrator.maximum_time_step = 1.0  # Limit max step size
    r.integrator.maximum_num_steps = 50000  # Increase max steps
    
    print("\n[1/3] Running environmental simulation (1500 minutes)...")
    print("  (Using stable integrator settings for stiff system)")
    
    try:
        result = r.simulate(0, 1500, 3000)
        print("✓ Simulation complete\n")
    except Exception as e:
        print(f"✗ Simulation failed: {e}")
        print("  Trying with even tighter tolerances...")
        r.integrator.absolute_tolerance = 1e-12
        r.integrator.maximum_num_steps = 100000
        result = r.simulate(0, 1500, 3000)
        print("✓ Simulation complete (second attempt)\n")
    
    return result, r


def plot_environmental_results(result):
    """Plot biosensor performance under environmental stress"""
    
    time = result[:, 0]
    
    # Biosensor variables
    sclerostin = result[:, 1]
    lrp6_free = result[:, 2]
    sclerostin_lrp6 = result[:, 3]
    wnt_active = result[:, 4]
    beta_catenin = result[:, 5]
    mrna = result[:, 6]
    reporter_inactive = result[:, 7]
    reporter_active = result[:, 8]
    
    # Environmental variables
    glucose = result[:, 9]
    atp = result[:, 10]
    amino_acids = result[:, 11]
    ros = result[:, 12]
    stress_protein = result[:, 13]
    temperature = result[:, 14]
    ph_level = result[:, 15]
    
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(4, 3, hspace=0.35, wspace=0.3)
    
    # Row 1: Environmental conditions
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(time, sclerostin, 'b-', linewidth=2.5)
    ax1.set_ylabel('Sclerostin (nM)', fontsize=11)
    ax1.set_title('INPUT: Patient Samples', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0.2, color='orange', linestyle='--', alpha=0.5, label='Clinical threshold')
    ax1.legend(fontsize=9)
    
    ax2 = fig.add_subplot(gs[0, 1])
    ax2_temp = ax2.twinx()
    ax2.plot(time, temperature, 'red', linewidth=2.5, label='Temp')
    ax2_temp.plot(time, ph_level, 'blue', linewidth=2.5, label='pH')
    ax2.set_ylabel('Temperature (°C)', fontsize=11, color='red')
    ax2_temp.set_ylabel('pH', fontsize=11, color='blue')
    ax2.set_title('ENVIRONMENT: Temp & pH', fontsize=12, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2_temp.tick_params(axis='y', labelcolor='blue')
    ax2.grid(True, alpha=0.3)
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_temp.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=9)
    
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(time, glucose, 'green', linewidth=2.5)
    ax3.set_ylabel('Glucose (nM)', fontsize=11)
    ax3.set_title('NUTRIENTS: Glucose', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Row 2: Cellular resources
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(time, atp, 'purple', linewidth=2.5)
    ax4.set_ylabel('ATP (nM)', fontsize=11)
    ax4.set_title('ENERGY: ATP Levels', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='Low energy')
    ax4.legend(fontsize=9)
    
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.plot(time, amino_acids, 'brown', linewidth=2.5)
    ax5.set_ylabel('Amino Acids (nM)', fontsize=11)
    ax5.set_title('BUILDING BLOCKS: AAs', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.plot(time, ros, 'red', linewidth=2.5, label='ROS')
    ax6.plot(time, stress_protein, 'orange', linewidth=2.5, label='Stress protein')
    ax6.set_ylabel('Concentration (nM)', fontsize=11)
    ax6.set_title('STRESS: Oxidative Damage', fontsize=12, fontweight='bold')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    
    # Row 3: Signaling cascade
    ax7 = fig.add_subplot(gs[2, 0])
    ax7.plot(time, lrp6_free, 'g-', linewidth=2, label='Free LRP6')
    ax7.plot(time, sclerostin_lrp6, 'r-', linewidth=2, label='Scl-LRP6')
    ax7.set_ylabel('Receptors (nM)', fontsize=11)
    ax7.set_title('DETECTION: Receptor States', fontsize=12, fontweight='bold')
    ax7.legend(fontsize=9)
    ax7.grid(True, alpha=0.3)
    
    ax8 = fig.add_subplot(gs[2, 1])
    ax8.plot(time, wnt_active, 'orange', linewidth=2.5, label='Wnt')
    ax8.plot(time, beta_catenin, 'purple', linewidth=2.5, label='β-catenin')
    ax8.set_ylabel('Concentration (nM)', fontsize=11)
    ax8.set_title('TRANSDUCTION: Signaling', fontsize=12, fontweight='bold')
    ax8.legend(fontsize=9)
    ax8.grid(True, alpha=0.3)
    
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.plot(time, mrna, 'brown', linewidth=2.5)
    ax9.set_ylabel('mRNA (nM)', fontsize=11)
    ax9.set_title('EXPRESSION: mRNA', fontsize=12, fontweight='bold')
    ax9.grid(True, alpha=0.3)
    
    # Row 4: Reporter output
    ax10 = fig.add_subplot(gs[3, 0])
    ax10.plot(time, reporter_inactive, 'gray', linewidth=2, label='Inactive')
    ax10.plot(time, reporter_active, 'lime', linewidth=2.5, label='Active')
    ax10.set_xlabel('Time (minutes)', fontsize=11)
    ax10.set_ylabel('Reporter (nM)', fontsize=11)
    ax10.set_title('OUTPUT: GFP Maturation', fontsize=12, fontweight='bold')
    ax10.legend(fontsize=9)
    ax10.grid(True, alpha=0.3)
    
    # Integrated view
    ax11 = fig.add_subplot(gs[3, 1:])
    
    # Normalize for overlay
    def safe_norm(x):
        max_x = np.max(np.abs(x))
        return x / max_x if max_x > 0 else x
    
    ax11.plot(time, safe_norm(sclerostin), 'b-', linewidth=2, label='Sclerostin (input)', alpha=0.8)
    ax11.plot(time, safe_norm(reporter_active), 'lime', linewidth=2.5, label='Reporter (output)', alpha=0.9)
    ax11.plot(time, safe_norm(atp - 50), 'purple', linewidth=1.5, label='ATP (energy)', alpha=0.6)
    ax11.plot(time, safe_norm(ros), 'red', linewidth=1.5, label='ROS (stress)', alpha=0.6)
    
    ax11.set_xlabel('Time (minutes)', fontsize=12)
    ax11.set_ylabel('Normalized Level', fontsize=12)
    ax11.set_title('INTEGRATED: Biosensor vs Environment', fontsize=13, fontweight='bold')
    ax11.legend(loc='upper right', fontsize=10)
    ax11.grid(True, alpha=0.3)
    
    # Event annotations
    events = [(100, 'Low\nScl'), (300, 'Cold'), (500, 'High\nScl'), 
              (700, 'pH'), (900, 'Starve'), (1100, 'ROS')]
    
    for t, label in events:
        ax11.axvline(x=t, color='gray', linestyle=':', alpha=0.3)
        ax11.text(t, 0.9, label, fontsize=8, ha='center', bbox=dict(boxstyle='round', 
                  facecolor='wheat', alpha=0.3))
    
    plt.suptitle('Biosensor Performance Under Environmental Stress',
                fontsize=15, fontweight='bold')
    plt.show()
    
    return fig


def analyze_robustness(result):
    """Quantify biosensor performance under stress"""
    
    time = result[:, 0]
    sclerostin = result[:, 1]
    reporter_active = result[:, 8]
    atp = result[:, 10]
    ros = result[:, 12]
    temperature = result[:, 14]
    
    print("\n" + "="*70)
    print(" ENVIRONMENTAL ROBUSTNESS ANALYSIS")
    print("="*70)
    
    conditions = {
        'Baseline (optimal)': (time >= 50) & (time <= 100),
        'Low sclerostin': (time >= 150) & (time <= 250),
        'Cold shock (30°C)': (time >= 330) & (time <= 380),
        'High sclerostin': (time >= 550) & (time <= 650),
        'pH stress (6.8)': (time >= 730) & (time <= 780),
        'Glucose starvation': (time >= 930) & (time <= 980),
        'Oxidative stress': (time >= 1130) & (time <= 1180)
    }
    
    print(f"\n{'Condition':<25} {'Reporter':<12} {'ATP':<10} {'ROS':<10} {'Status'}")
    print("-" * 70)
    
    for condition, mask in conditions.items():
        rep = np.mean(reporter_active[mask])
        atp_avg = np.mean(atp[mask])
        ros_avg = np.mean(ros[mask])
        
        # Performance score
        atp_ok = atp_avg > 2000
        ros_ok = ros_avg < 2.0
        rep_ok = rep > 0.1 if np.mean(sclerostin[mask]) < 0.3 else rep < 5.0
        
        status = "✓ Good" if (atp_ok and ros_ok) else "⚠ Stress" if atp_ok else "✗ Fail"
        
        print(f"{condition:<25} {rep:>8.2f} nM  {atp_avg:>6.1f} nM  {ros_avg:>6.2f} nM  {status}")
    
    print("="*70 + "\n")


def parameter_sensitivity_analysis(model_obj, baseline_result):
    """
    Quick sensitivity analysis: vary key parameters ±20%
    """
    
    print("\n" + "="*70)
    print(" PARAMETER SENSITIVITY ANALYSIS")
    print("="*70)
    
    # Key parameters to test
    params_to_test = [
        ('k_on_base', 0.5),
        ('k_wnt_activation_base', 0.05),
        ('k_transcription_base', 0.08),
        ('k_atp_production', 2.0),
        ('k_ros_scavenge', 0.08)
    ]
    
    print(f"\n{'Parameter':<25} {'Baseline':<12} {'-20%':<12} {'+20%':<12} {'Sensitivity'}")
    print("-" * 70)
    
    for param_name, baseline_value in params_to_test:
        # Get baseline reporter output
        baseline_reporter = np.mean(baseline_result[200:300, 8])  # Low sclerostin period
        
        # Test -20%
        model_obj.reset()
        model_obj[param_name] = baseline_value * 0.8
        result_low = model_obj.simulate(0, 1500, 1000)
        reporter_low = np.mean(result_low[200:300, 8])
        
        # Test +20%
        model_obj.reset()
        model_obj[param_name] = baseline_value * 1.2
        result_high = model_obj.simulate(0, 1500, 1000)
        reporter_high = np.mean(result_high[200:300, 8])
        
        # Calculate sensitivity
        sensitivity = (reporter_high - reporter_low) / (0.4 * baseline_value * baseline_reporter)
        sensitivity_pct = abs(sensitivity * 100)
        
        # Reset to baseline
        model_obj.reset()
        model_obj[param_name] = baseline_value
        
        status = "HIGH" if sensitivity_pct > 50 else "MED" if sensitivity_pct > 20 else "LOW"
        
        print(f"{param_name:<25} {baseline_reporter:>8.2f} nM  {reporter_low:>8.2f} nM  "
              f"{reporter_high:>8.2f} nM  {status} ({sensitivity_pct:.1f}%)")
    
    print("="*70 + "\n")



def main():
    """Run complete environmental simulation"""
    
    result, model = run_environmental_simulation()
    
    print("[2/4] Generating environmental analysis plots...")
    plot_environmental_results(result)
    
    print("[3/4] Analyzing biosensor robustness...")
    analyze_robustness(result)
    
    print("[4/4] Running parameter sensitivity analysis...")
    parameter_sensitivity_analysis(model, result)
    
    print("✓ Environmental simulation complete!\n")
    
    return result, model


if __name__ == "__main__":
    result, model = main()
