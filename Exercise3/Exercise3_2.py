import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── Konstanten ─────────────────────────────────────────────────────────────────
# GM = 4π²: aus Keplers 3. Gesetz T²=4π²a³/(GM)
# mit T=1 YR, a=1 AU (Erde) → GM = 4π² AU³/YR²
GM      = 4 * np.pi**2
EPSILON = 1e-3   # max. relativer Fehler ε pro Schritt
S1      = 0.9    # Sicherheitsfaktor 1 (Gl. 2.61, Vorlesung)
S2      = 4.0    # Sicherheitsfaktor 2 (Gl. 2.61, Vorlesung)

# ── Bewegungsgleichungen ───────────────────────────────────────────────────────
def f(state):
    """d/dt [x, y, vx, vy] = [vx, vy, ax, ay]"""
    x, y, vx, vy = state
    r3 = (x**2 + y**2)**1.5
    return np.array([vx, vy, -GM*x/r3, -GM*y/r3])

# ── Einzelner RK4-Schritt ──────────────────────────────────────────────────────
def rk4_step(state, tau):
    k1 = f(state)
    k2 = f(state + 0.5*tau*k1)
    k3 = f(state + 0.5*tau*k2)
    k4 = f(state +     tau*k3)
    return state + (tau/6.0) * (k1 + 2*k2 + 2*k3 + k4)

# ── Adaptiver RK4-Schritt (genau nach Vorlesung, Abschnitt 2.8) ───────────────
def adaptive_step(state, tau):
    """
    Vorlesung Gl. 2.60 / 2.61 / 2.62:
      x̄  = ein voller Schritt der Größe τ
      x̃  = zwei halbe Schritte der Größe τ/2
      Fehlermaß (multidimensional, Gl. 2.62):
        δ = max_i { |x̄_i - x̃_i| / ( (|x̄_i|+|x̃_i|)/2 + Δ ) }
      Akzeptiert wird x̃, KEIN Richardson-verbesserter Wert.
    """
    x_bar  = rk4_step(state, tau)                          # ein voller Schritt  (x̄)
    x_tilde = rk4_step(rk4_step(state, tau/2), tau/2)      # zwei halbe Schritte (x̃)

    Delta  = 1e-12   # kleine Zahl gegen Division durch 0 (Gl. 2.62)
    denom  = 0.5 * (np.abs(x_bar) + np.abs(x_tilde)) + Delta   # (|x̄|+|x̃|)/2 + Δ
    delta  = np.max(np.abs(x_bar - x_tilde) / denom)            # Fehlermaß δ (Gl. 2.62)

    return x_tilde, delta   # x̃ wird akzeptiert (Vorlesung: KEIN /15, KEIN Richardson)

# ── Schrittweitenanpassung (genau nach Vorlesung, Gl. 2.61) ───────────────────
def update_tau(tau, delta):
    """
    τ' = τ * (ε/δ)^(1/5)        [Gl. 2.61]
    τ_new = S1 * τ'
    Grenzen: τ/S2 ≤ τ_new ≤ S2*τ  [Vorlesung, S1=0.9, S2=4]
    """
    if delta == 0:
        return tau * S2
    tau_prime = tau * (EPSILON / delta)**0.2
    tau_new   = S1 * tau_prime
    if   tau_new > S2 * tau:   return S2 * tau     # obere Grenze
    elif tau_new < tau / S2:   return tau / S2      # untere Grenze (= τ/4 !)
    else:                       return tau_new

# ── Simulation ─────────────────────────────────────────────────────────────────
def simulate(t_max, tau0=0.1):
    state = np.array([1.0, 0.0, 0.0, np.pi/2])   # Aphelion: r=1AU, v=π/2 AU/YR
    tau   = tau0
    t     = 0.0

    states, times, taus, radii = [state.copy()], [0.0], [], []
    max_reductions = 50   # Warnung nach zu vielen Verkleinerungen

    while t < t_max:
        attempts = 0
        while True:
            x_new, delta = adaptive_step(state, tau)

            if delta <= EPSILON or tau < 1e-9:
                # Schritt akzeptiert
                state  = x_new
                t     += tau
                states.append(state.copy())
                times.append(t)
                r = np.sqrt(state[0]**2 + state[1]**2)
                taus.append(tau)
                radii.append(r)
                tau = update_tau(tau, delta)   # Schrittweite für nächsten Schritt
                break
            else:
                # Schritt abgelehnt → τ verkleinern
                tau = update_tau(tau, delta)
                attempts += 1
                if attempts > max_reductions:
                    print(f"Warnung: mehr als {max_reductions} Verkleinerungen bei t={t:.4f}")
                    break

    return np.array(states), np.array(times), np.array(taus), np.array(radii)

# ── Analytische Lösung ─────────────────────────────────────────────────────────
# Drehimpuls: L = r_a * v_a = 1 * π/2 = π/2
# Energie:    E = v²/2 - GM/r = π²/8 - 4π²
# Umkehrpunkte: L²/(2r²) - GM/r = E  →  31r² - 32r + 1 = 0
# → r_Aphel = 1, r_Perihel = 1/31
r_aph_exact  = 1.0
r_per_exact  = 1.0 / 31.0
e_exact      = (r_aph_exact - r_per_exact) / (r_aph_exact + r_per_exact)  # = 15/16
a_exact      = (r_aph_exact + r_per_exact) / 2.0
T_exact      = a_exact**1.5   # Kepler: T² = a³ → T = a^(3/2) in diesen Einheiten

print("=" * 52)
print("  ANALYTISCHE LÖSUNG")
print("=" * 52)
print(f"  Aphel-Abstand   r_a = {r_aph_exact:.6f} AU")
print(f"  Perihel-Abstand r_p = {r_per_exact:.6f} AU  (= 1/31)")
print(f"  Halbachse        a  = {a_exact:.6f} AU")
print(f"  Exzentrizität    e  = {e_exact:.6f}  (= 15/16)")
print(f"  Umlaufzeit       T  = {T_exact:.6f} YR")
print()

# ── Simulation ausführen (ca. 2 Perioden für guten log-log Plot) ───────────────
states, times, taus, radii = simulate(t_max=2*T_exact, tau0=0.1)

x, y   = states[:, 0], states[:, 1]
r      = np.sqrt(x**2 + y**2)
idx_per = np.argmin(r)

r_aph_num = 1.0   # exakt durch Anfangsbedingung
r_per_num = r[idx_per]
e_num     = (r_aph_num - r_per_num) / (r_aph_num + r_per_num)

print("=" * 52)
print("  NUMERISCHE ERGEBNISSE")
print("=" * 52)
print(f"  Aphel-Abstand   r_a = {r_aph_num:.6f} AU  (Anfangsbedingung)")
print(f"  Perihel-Abstand r_p = {r_per_num:.6f} AU")
print(f"  Exzentrizität    e  = {e_num:.6f}")
print(f"  Anzahl Schritte     = {len(taus)}")
print()

# ── Energieerhaltung ───────────────────────────────────────────────────────────
def energy(s):
    x, y, vx, vy = s
    return 0.5*(vx**2+vy**2) - GM/np.sqrt(x**2+y**2)
E0, E1 = energy(states[0]), energy(states[-1])
print(f"  Energieerhaltung: rel. Abweichung = {abs(E1-E0)/abs(E0)*100:.5f} %")
print()

# ── Log-Log Fit: τ vs r (Keplers 3. Gesetz → Steigung ≈ 3/2) ─────────────────
log_r   = np.log10(radii)
log_tau = np.log10(taus)
slope, intercept = np.polyfit(log_r, log_tau, 1)
print(f"  Log-Log-Steigung τ vs r : {slope:.4f}")
print(f"  Erwartung (Kepler 3)    : 1.5  (da τ ∝ r^(3/2))")
print("=" * 52)

# ── Plots ──────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10))
fig.patch.set_facecolor('#0d1117')
gs  = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.35)
ax_orbit  = fig.add_subplot(gs[:, 0])
ax_tau    = fig.add_subplot(gs[0, 1])
ax_loglog = fig.add_subplot(gs[1, 1])
lkw = dict(color='#c9d1d9', fontsize=11)

# — Bahn —
ax_orbit.set_facecolor('#0d1117')
sc = ax_orbit.scatter(x, y, c=r, cmap='plasma', s=5, zorder=3)
cb = plt.colorbar(sc, ax=ax_orbit, pad=0.02)
cb.set_label('Abstand [AU]', color='#c9d1d9', fontsize=10)
plt.setp(cb.ax.yaxis.get_ticklabels(), color='#8b949e')
ax_orbit.scatter(0, 0, color='#fbbf24', s=250, zorder=5, label='Sonne', marker='*')
ax_orbit.scatter(x[0], y[0], color='#60a5fa', s=90, zorder=6,
                 label=f'Aphel  r = {r_aph_num:.4f} AU')
ax_orbit.scatter(x[idx_per], y[idx_per], color='#f87171', s=90, zorder=6,
                 label=f'Perihel r = {r_per_num:.5f} AU')
ax_orbit.set_xlabel('x [AU]', **lkw); ax_orbit.set_ylabel('y [AU]', **lkw)
ax_orbit.set_title(f'Kometenbahn  (e = {e_num:.5f},  e_exakt = {e_exact:.5f})',
                   color='#c9d1d9', fontsize=12, fontweight='bold')
ax_orbit.set_aspect('equal')
ax_orbit.tick_params(colors='#8b949e')
ax_orbit.legend(fontsize=9, facecolor='#161b22', labelcolor='#c9d1d9', edgecolor='#30363d')
for sp in ax_orbit.spines.values(): sp.set_color('#30363d')

# — τ vs Zeit —
ax_tau.set_facecolor('#0d1117')
ax_tau.plot(times[1:], taus, color='#34d399', lw=0.9, alpha=0.9)
ax_tau.set_xlabel('Zeit [YR]', **lkw); ax_tau.set_ylabel('Schrittweite τ [YR]', **lkw)
ax_tau.set_title('Adaptive Schrittweite τ(t)', color='#c9d1d9', fontsize=11, fontweight='bold')
ax_tau.tick_params(colors='#8b949e')
for sp in ax_tau.spines.values(): sp.set_color('#30363d')

# — Log-Log: τ vs r —
ax_loglog.set_facecolor('#0d1117')
ax_loglog.scatter(radii, taus, s=12, alpha=0.55, color='#818cf8', label='Simulationsdaten')
r_fit   = np.logspace(np.log10(radii.min()), np.log10(radii.max()), 300)
tau_fit = 10**intercept * r_fit**slope
ax_loglog.loglog(r_fit, tau_fit, color='#f97316', lw=2.5,
                 label=f'Fit: Steigung = {slope:.3f}\n(Kepler: 1.5)')
ax_loglog.set_xlabel('Abstand r [AU]', **lkw); ax_loglog.set_ylabel('Schrittweite τ [YR]', **lkw)
ax_loglog.set_title('Log-Log: τ vs r', color='#c9d1d9', fontsize=11, fontweight='bold')
ax_loglog.legend(fontsize=9, facecolor='#161b22', labelcolor='#c9d1d9', edgecolor='#30363d')
ax_loglog.tick_params(colors='#8b949e')
for sp in ax_loglog.spines.values(): sp.set_color('#30363d')


plt.savefig('comet_orbit_new.png', dpi=150,
            bbox_inches='tight', facecolor='#0d1117')
print("\nFigur gespeichert.")