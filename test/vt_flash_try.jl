using CubicEoS

mix = load(BrusilovskyEoSMixture;
    # names=("methane",),
    # names=("methane", "n-pentane"),
    # names=("nitrogen", "methane", "propane", "n-decane"),
    # names=("nitrogen","carbon_dioxide","hydrogen_sulfide","methane","ethane","propane","n-butane","n-pentane"),
    # names=("nitrogen","carbon_dioxide","hydrogen_sulfide","methane","ethane","propane","i-butane","n-butane","i-pentane","n-pentane","n-hexane","n-heptane","n-octane","n-nonane","n-decane","undecane","dodecane","tridecane","tetradecane","pentadecane","hexadecane","heptadecane","octadecane","nonadecane","icosane"),
    names=("carbon_dioxide", "n-decane"),
)

# χ = [1.0]
χ = [0.547413, 0.452587]
# χ = [0.2463, 0.2208, 0.2208, 0.3121]
# χ = [0.00325,0.01556,0.03329,0.82829,0.08552,0.01725,0.01015,0.00669]
# χ = [0.00325,0.01556,0.03329,0.82829,0.05850,0.02702,0.00543,0.00562,0.00248,0.00198,0.00174,0.00179,0.00275,0.00233,0.00191,0.00137,0.00112,0.00098,0.00087,0.00076,0.00053,0.00042,0.00039,0.00031,0.00131]


V = 1e-6
N = (4000 * V) .* χ  # 8500 phase transition c1c5, c1 4000
RT = CubicEoS.GAS_CONSTANT_SI * 250  # 170 c1, 371 c1c5

state = CubicEoS.vt_flash(mix, N, V, RT)
dump(state)
