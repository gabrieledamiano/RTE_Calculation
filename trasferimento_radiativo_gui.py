#!/usr/bin/env python3
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

# =============================================================================
# FUNZIONI UTILI
# =============================================================================

def blackbodyn(T, numerionda):
    """
    Calcola la radianza del corpo nero per una data temperatura T (K)
    e un vettore di numeri d'onda (cm⁻¹).
    """
    c_local = 2.998e8       # m/s
    h_local = 6.626e-34     # J*s
    k_local = 1.381e-23     # J/K
    c1 = 2 * h_local * c_local**2 * 1e8
    c2 = h_local * c_local / k_local * 1e2
    BB = c1 * numerionda**3 * np.exp(-c2 * numerionda / T) / (1.0 - np.exp(-c2 * numerionda / T))
    return BB

def lorentzian(w0, S, a0, w, dwm=0, T=296, p=1013.25):
    """
    Calcola il profilo lorentziano per un insieme di righe spettrali.

    Parametri:
      w0  : array dei numeri d'onda centrali (cm⁻¹)
      S   : array dell'intensità (da dividere per π)
      a0  : array della larghezza a metà altezza (cm⁻¹)
      w   : array dei numeri d'onda dove calcolare il profilo (cm⁻¹)
      dwm : eventuale taglio (default 0, nessun taglio)
      T   : temperatura (K)
      p   : pressione (hPa)
    
    Restituisce:
      f   : profilo di assorbimento (array 1D, lunghezza pari a len(w))
    """
    w0 = np.array(w0).reshape(-1, 1)   # (n_linee, 1)
    S = np.array(S).reshape(-1, 1) / np.pi
    a0 = np.array(a0).reshape(-1, 1)
    w = np.array(w).reshape(1, -1)       # (1, n_w)
    
    T0 = 296.0
    p0 = 1013.25
    a = (T0 / T)**0.5 * (p / p0) * a0
    
    f = np.zeros((1, w.shape[1]))
    chunk_size = 1000
    n_lines = w0.shape[0]
    num_chunks = int(np.ceil(n_lines / chunk_size))
    
    for i in range(num_chunks):
        idx_start = i * chunk_size
        idx_end = min((i+1)*chunk_size, n_lines)
        w0_chunk = w0[idx_start:idx_end]
        S_chunk = S[idx_start:idx_end]
        a_chunk = a0[idx_start:idx_end]
        
        DW = w - w0_chunk 
        if dwm > 0:
            mask = np.abs(DW) <= dwm
        else:
            mask = np.ones_like(DW, dtype=bool)
            
        f_chunk = (S_chunk * a_chunk) / (DW**2 + a_chunk**2) * mask
        f += np.sum(f_chunk, axis=0, keepdims=True)
    
    return f.flatten()

def calcola_altezze_strati(p, t):
    """
    Calcola l'altezza degli strati atmosferici (in Km) a partire dai vettori
    delle pressioni (p in hPa) e delle temperature (t in K).
    """
    g0 = 9.80665
    Na_local = 6.022142e23
    Re_local = 6371       # Km
    k_local = 1.38065e-23
    mair_local = 28.964
    R_const = k_local * Na_local
    p = np.array(p).flatten()
    t = np.array(t).flatten()
    
    Pu = p[1:] * 100.0   # conversione in Pa
    Tu = t[1:]
    Pl = p[:-1] * 100.0  # conversione in Pa
    Tl = t[:-1]
    
    Z = np.zeros(len(p))
    Tm = np.zeros(len(p)-1)
    Z[0] = 0.0
    for i in range(1, len(p)):
        tup = Tu[i-1]
        tlo = Tl[i-1]
        pup = Pu[i-1]
        plo = Pl[i-1]
        Tm[i-1] = tup + (tlo - tup)/np.log(plo/pup) * ((plo*(np.log(plo)-1) - pup*(np.log(pup)-1))/(plo-pup) - np.log(pup))
        G = g0 * Re_local**2 / (Re_local + Z[i-1])**2
        Z[i] = Z[i-1] - np.log(pup/plo) * R_const * Tm[i-1] / (mair_local * G)
    return Z

def calcola_medie(p, t):
    """
    Calcola la pressione e la temperatura medie per ciascun strato.
    
    Input:
      p : array delle pressioni (hPa) (lunghezza N+1)
      t : array delle temperature (K) (lunghezza N+1)
      
    Restituisce:
      pm : array delle pressioni medie (hPa) (lunghezza N)
      tm : array delle temperature medie (K) (lunghezza N)
    """
    p = np.array(p).flatten()
    t = np.array(t).flatten()
    n = len(p) - 1
    pm = np.zeros(n)
    tm = np.zeros(n)
    for i in range(n):
        pm[i] = (p[i+1] + p[i]) / 2.0
        tm[i] = (t[i+1] + t[i]) / 2.0
    return pm, tm

def calcola_gas_col(t, p, gas_mr, Z, Mgas=None):
    """
    Calcola la densità colonnare di un gas (in molecole/m²) in atmosfera.
    
    Parametri:
      t      : array delle temperature (K)
      p      : array delle pressioni (hPa)
      gas_mr : array del mixing ratio (se Mgas è specificato, si assume in g/kg, altrimenti frazione molare)
      Z      : array delle altezze (Km) (da calcola_altezze_strati)
      Mgas   : massa molecolare del gas (g/mol) (opzionale)
      
    Restituisce:
      gas_col : array della densità colonnare per ciascun strato.
    """
    pm, tm = calcola_medie(p, t)
    pm = pm * 100.0  # conversione da hPa a Pa
    R_const = 8.314
    Na_local = 6.022e23
    n_levels = len(p)
    gas_col = np.zeros(n_levels - 1)
    gas_mr = np.array(gas_mr).flatten()
    if Mgas is not None:
        # conversione da g/kg a frazione molare
        gas_kg = gas_mr * 0.001
        M_air = 28.97
        epsilon = Mgas / M_air
        chi = gas_kg / (gas_kg + epsilon)
    else:
        chi = gas_mr
    for i in range(n_levels - 1):
        dz = (Z[i+1] - Z[i])
        # Converte dz da Km a m (×1000) e moltiplica per 1e-4 (come nel codice MATLAB)
        moli_aria_per_m2 = (pm[i] / (R_const * tm[i])) * dz * 1000 * 1e-4
        gas_col[i] = chi[i] * moli_aria_per_m2 * Na_local
    return gas_col

# =============================================================================
# GUI PER IL CALCOLO DELLA RADIANZA
# =============================================================================

class RadianceGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Calcolo Radianza Atmosferica")
        self.geometry("1100x750")
        self.dark_mode = True  # Variabile per tenere traccia del tema
        self.configure(bg='#1C2331')  # Deep navy blue, reminiscent of night sky
        
        # Colori tema con nuances di celeste chiaro e design moderno
        self.colors = {
            'dark': {
                'bg': '#1C2331',  # Deep navy blue
                'fg': '#87CEEB',  # Light sky blue
                'accent_primary': '#87CEEB',  # Light sky blue
                'accent_secondary': '#4682B4',  # Steel blue
                'button_bg': '#34495E',  # Muted blue-gray
                'button_fg': '#E6F2FF',  # Very light blue
                'button_hover': '#5F9EA0',  # Cadet blue
                'text_primary': '#87CEEB',  # Light sky blue
                'text_secondary': '#B0E0E6',  # Powder blue
                'frame_bg': '#1C2331',
                'border_color': '#2C3E50',
                'shadow_color': 'rgba(0,0,0,0.3)'
            },
            'light': {
                'bg': '#F0F8FF',  # Alice Blue
                'fg': '#4682B4',  # Steel blue
                'accent_primary': '#4682B4',  # Steel blue
                'accent_secondary': '#87CEEB',  # Light sky blue
                'button_bg': '#2980B9',  # Bright blue
                'button_fg': '#FFFFFF',  # White
                'button_hover': '#3498DB',  # Lighter blue
                'text_primary': '#2C3E50',  # Dark slate
                'text_secondary': '#34495E',  # Muted blue-gray
                'frame_bg': '#F0F8FF',
                'border_color': '#BDC3C7',
                'shadow_color': 'rgba(0,0,0,0.1)'
            }
        }

        # Stile più sofisticato con ombreggiature e bordi arrotondati
        self.style = tk.ttk.Style()
        self.style.theme_use('clam')
        
        # Configurazione stile pulsanti
        self.style.configure('TButton', 
                             font=('Arial', 10, 'bold'),
                             borderwidth=0,
                             relief='flat')
        self.style.map('TButton',
                       background=[('active', self.get_current_color('button_hover'))])

        # Frame principale con ombra e bordi arrotondati
        self.main_frame = tk.Frame(self, 
                                   bg=self.get_current_color('bg'),
                                   bd=0,
                                   highlightthickness=0)
        self.main_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=20)

        # Pulsante tema con design moderno
        self.theme_button = tk.Button(
            self.main_frame, 
            text="☀️" if self.dark_mode else "🌙", 
            font=('Segoe UI Emoji', 16),
            command=self.toggle_theme,
            borderwidth=0,
            highlightthickness=0,
            padx=10,
            pady=5,
            bg=self.get_current_color('bg'),
            activebackground=self.get_current_color('bg'),
            fg=self.get_current_color('accent_primary'),
            activeforeground=self.get_current_color('accent_secondary')
        )
        self.theme_button.pack(side=tk.TOP, anchor=tk.E, padx=10, pady=5)

        # Intestazione con design più elegante
        header_frame = tk.Frame(self.main_frame, 
                                bg=self.get_current_color('bg'),
                                bd=0)
        header_frame.pack(fill=tk.X, pady=(10, 20))
        
        # Icone atmosferiche con animazione subtile
        header_icon = tk.Label(header_frame, 
                               text="🌍 ☀️ ☁️ 🌈", 
                               font=('Segoe UI Emoji', 40), 
                               bg=self.get_current_color('bg'), 
                               fg=self.get_current_color('accent_primary'))
        header_icon.pack(side=tk.TOP, anchor='center', pady=(0, 10))
        
        # Titolo con design più moderno
        header_label = tk.Label(header_frame, 
                                text="Calcolo Radianza Atmosferica", 
                                font=('Arial', 18, 'bold'), 
                                bg=self.get_current_color('bg'), 
                                fg=self.get_current_color('text_primary'),
                                anchor='center',
                                justify='center')
        header_label.pack(side=tk.TOP, expand=True, pady=(0, 5))
        
        # Descrizione scientifica con stile più elegante
        description_label = tk.Label(self.main_frame, 
                                     text="Simulazione dell'emissione radiativa\nnell'atmosfera terrestre", 
                                     font=('Arial', 11, 'italic'), 
                                     bg=self.get_current_color('bg'), 
                                     fg=self.get_current_color('text_secondary'),
                                     anchor='center',  
                                     justify='center')
        description_label.pack(pady=(0, 20))

        # Frame controlli con design più moderno
        control_frame = tk.Frame(self.main_frame, 
                                 bg=self.get_current_color('bg'), 
                                 bd=1, 
                                 relief=tk.FLAT)
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=20, pady=10)

        # Stile pulsanti più accattivante
        self.load_button = tk.ttk.Button(control_frame, 
                                         text="Carica Atmosfera", 
                                         command=self.load_atmosphere,
                                         style='Accent.TButton')
        self.load_button.pack(side=tk.LEFT, padx=10)
        
        self.calc_button = tk.ttk.Button(control_frame, 
                                         text="Calcola Radianza", 
                                         command=self.calculate_radiance, 
                                         state=tk.DISABLED,
                                         style='Accent.TButton')
        self.calc_button.pack(side=tk.LEFT, padx=10)
        
        # Etichetta file con design più pulito
        self.file_label = tk.Label(control_frame, 
                                   text="Nessun file caricato", 
                                   font=('Arial', 10, 'italic'), 
                                   bg=self.get_current_color('bg'), 
                                   fg=self.get_current_color('text_secondary'))
        self.file_label.pack(side=tk.LEFT, padx=10)

        # Frame per l'area grafica con bordo che ricorda gli strati atmosferici
        self.canvas_frame = tk.Frame(self.main_frame, 
                                     bg=self.get_current_color('bg'), 
                                     bd=2, 
                                     relief=tk.GROOVE)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.atmosphere_data = None  # per memorizzare i dati caricati

    def update_theme_style(self):
        """Aggiorna lo stile in base al tema corrente"""
        current_theme = 'dark' if self.dark_mode else 'light'
        colors = self.colors[current_theme]
        
        # Aggiorna i colori della finestra e dei frame
        self.configure(bg=colors['bg'])
        if hasattr(self, 'main_frame'):
            self.main_frame.configure(bg=colors['frame_bg'])
            
            # Aggiorna il pulsante del tema
            if hasattr(self, 'theme_button'):
                self.theme_button.configure(
                    bg=colors['frame_bg'],
                    activebackground=colors['frame_bg'],
                    fg=colors['accent_primary']
                )
            
            # Aggiorna ricorsivamente tutti i widget
            def update_widget_colors(widget):
                if isinstance(widget, tk.Frame):
                    widget.configure(bg=colors['frame_bg'])
                elif isinstance(widget, tk.Label):
                    widget.configure(bg=colors['frame_bg'], fg=colors['text_secondary'])
                # Aggiorna ricorsivamente tutti i widget figli
                for child in widget.winfo_children():
                    update_widget_colors(child)
            
            update_widget_colors(self.main_frame)

    def get_current_color(self, key):
        """Ottiene il colore corrente in base al tema"""
        current_theme = 'dark' if self.dark_mode else 'light'
        return self.colors[current_theme][key]
    
    def toggle_theme(self):
        """Cambia il tema tra chiaro e scuro"""
        self.dark_mode = not self.dark_mode
        # Aggiorna l'icona del pulsante tema
        self.theme_button.configure(
            text="☀️" if self.dark_mode else "🌙",
            bg=self.get_current_color('bg'),
            activebackground=self.get_current_color('bg'),
            fg=self.get_current_color('accent_primary'),
            activeforeground=self.get_current_color('accent_secondary')
        )
        self.update_theme_style()
    
    def load_atmosphere(self):
        file_path = filedialog.askopenfilename(title="Seleziona file atmosfera (MAT)",
                                            filetypes=[("MAT files", "*.mat")])
        if file_path:
            try:
                self.atmosphere_data = loadmat(file_path)
                self.file_label.config(text=file_path)
                self.calc_button.config(state=tk.NORMAL)
                messagebox.showinfo("Successo", "File atmosfera caricato correttamente!")
                
                # Plot pressure vs temperature and mixing ratio side by side
                p = self.atmosphere_data.get('p').squeeze()
                t = self.atmosphere_data.get('t').squeeze()
                pm, tm = calcola_medie(p, t)
                
                figure = Figure(figsize=(12, 5))
                
                # Pressure vs Temperature subplot
                ax1 = figure.add_subplot(121)
                ax1.plot(tm, pm)
                ax1.set_xlabel('Temperatura Media (K)')
                ax1.set_ylabel('Pressione Media (hPa)')
                ax1.set_title('Pressione vs Temperatura')
                ax1.set_yscale('log')
                ax1.invert_yaxis()
                ax1.grid(True, which="both", ls="-", alpha=0.2)
                
                # Mixing Ratio vs Pressure subplot con 12 colori distinti
                ax2 = figure.add_subplot(122)
                mr_keys = [key for key in self.atmosphere_data.keys() if key.endswith('_mr') and not key.startswith('__')]
                
                # Importiamo i moduli necessari per il colormap
                import matplotlib.cm as cm
                import numpy as np
                # Creiamo una lista di 12 colori distinti usando il colormap 'tab20'
                colors = [cm.get_cmap('tab20')(i) for i in np.linspace(0, 1, 12)]
                
                # Assumiamo che ci siano 12 gas; in caso contrario, si può usare modulo o gestire diversamente
                for i, key in enumerate(mr_keys):
                    gas_mr = self.atmosphere_data[key].squeeze()
                    # Usa il colore corrispondente dalla lista
                    ax2.plot(gas_mr, pm, label=key.replace('_mr', ''), color=colors[i])
                
                ax2.set_xlabel('Mixing Ratio')
                ax2.set_ylabel('Pressione Media (hPa)')
                ax2.set_title('Mixing Ratio vs Pressione')
                ax2.set_xscale('log')
                ax2.set_yscale('log')
                ax2.invert_yaxis()
                ax2.legend()
                ax2.grid(True, which="both", ls="-", alpha=0.2)
                
                figure.tight_layout()
                
                # Pulisce il canvas precedente e inserisce il nuovo grafico
                for widget in self.canvas_frame.winfo_children():
                    widget.destroy()
                
                canvas = FigureCanvasTkAgg(figure, master=self.canvas_frame)
                canvas.draw()
                canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
                
            except Exception as e:
                messagebox.showerror("Errore", f"Impossibile caricare il file: {e}")


    def calculate_radiance(self):
        if self.atmosphere_data is None:
            messagebox.showerror("Errore", "Caricare prima il file atmosfera!")
            return
        try:
            # Crea una finestra di progresso
            progress_window = tk.Toplevel(self)
            progress_window.title("Calcolo in corso...")
            progress_window.geometry("300x150")
            progress_window.transient(self)  # Rende la finestra dipendente dalla principale
            
            # Centra la finestra
            progress_window.geometry("+%d+%d" % (
                self.winfo_x() + self.winfo_width()/2 - 150,
                self.winfo_y() + self.winfo_height()/2 - 75))
            
            # Label per il progresso
            progress_label = tk.Label(progress_window, text="Calcolo della radianza in corso...", pady=10)
            progress_label.pack()
            
            # Barra di progresso
            progress_bar = ttk.Progressbar(progress_window, mode='determinate', length=200)
            progress_bar.pack(pady=20)
            
            def update_progress(value, text=""):
                progress_bar['value'] = value
                if text:
                    progress_label.config(text=text)
                progress_window.update()
            
            # Estrae le variabili p e t dal file MAT
            update_progress(10, "Estrazione dati dall'atmosfera...")
            p = self.atmosphere_data.get('p')
            t = self.atmosphere_data.get('t')
            if p is None or t is None:
                progress_window.destroy()
                messagebox.showerror("Errore", "Le variabili 'p' e 't' non sono presenti nel file!")
                return
            p = p.squeeze()
            t = t.squeeze()
            
            # Definisce l'intervallo dei numeri d'onda (cm⁻¹)
            update_progress(20, "Definizione intervallo spettrale...")
            w_min = 900
            w_max = 1100
            w = np.arange(w_min, w_max + 0.01, 0.01)
            
            # Calcola le altezze degli strati
            update_progress(30, "Calcolo altezze degli strati...")
            Z = calcola_altezze_strati(p, t)
            
            # Calcola la densità colonnare dei gas
            update_progress(40, "Calcolo densità colonnare dei gas...")
            mr_keys = [key for key in self.atmosphere_data.keys() if key.endswith('_mr') and not key.startswith('__')]
            gas_columns = {}
            for key in mr_keys:
                gas_mr = self.atmosphere_data[key].squeeze()
                new_key = key.replace('_mr', '_col')
                if key == 'H2O_mr':
                    gas_columns[new_key] = calcola_gas_col(t, p, gas_mr, Z, Mgas=18.01)
                else:
                    gas_columns[new_key] = calcola_gas_col(t, p, gas_mr, Z)
            
            # Estrae H2O_col
            H2O_col = gas_columns.get('H2O_col')
            
            # Carica i dati spettrali
            update_progress(50, "Caricamento dati spettrali...")
            h2o_data  = np.loadtxt('H2O_hitran.txt')
            co2_data  = np.loadtxt('CO2_hitran.txt')
            o3_data   = np.loadtxt('O3_hitran.txt')
            nh3_data  = np.loadtxt('NH3_hitran.txt')
            hno3_data = np.loadtxt('HNO3_hitran.txt')
            
            update_progress(60, "Calcolo parametri atmosferici...")
            pm, tm = calcola_medie(p, t)
            n_layers = len(pm)
            n_w = len(w)
            
            # Parametri angolari e di superficie
            solar_zenith_angle = 45 * np.pi / 180
            viewing_angle = 0 * np.pi / 180
            mu_s = np.cos(solar_zenith_angle)
            mu = np.cos(viewing_angle)
            surface_emissivity = 0.7
            surface_temperature = 300  # K
            
            # Inizializzazione delle matrici delle trasmittanze
            tau = np.ones((n_w, n_layers + 1))
            tau_mu = np.ones((n_w, n_layers + 1))

            update_progress(70, "Calcolo spessore ottico...")
            for l in range(n_layers - 1, -1, -1):
                TT = tm[l]
                PP = pm[l]
                H2O_spec = lorentzian(h2o_data[:, 2], h2o_data[:, 3], h2o_data[:, 4], w, 0, TT, PP)
                CO2_spec = lorentzian(co2_data[:, 2], co2_data[:, 3], co2_data[:, 4], w, 0, TT, PP)
                O3_spec  = lorentzian(o3_data[:, 2], o3_data[:, 3], o3_data[:, 4], w, 0, TT, PP)
                NH3_spec = lorentzian(nh3_data[:, 2], nh3_data[:, 3], nh3_data[:, 4], w, 0, TT, PP)
                HNO3_spec= lorentzian(hno3_data[:, 2], hno3_data[:, 3], hno3_data[:, 4], w, 0, TT, PP)
                
                opd = (H2O_spec * H2O_col[l] +
                       CO2_spec * gas_columns['CO2_col'][l] +
                       O3_spec  * gas_columns['O3_col'][l] +
                       NH3_spec * gas_columns['NH3_col'][l] +
                       HNO3_spec* gas_columns['HNO3_col'][l])
                
                tau[:, l] = tau[:, l+1] * np.exp(-opd)
                tau_mu[:, l] = tau_mu[:, l+1] * np.exp(-opd / mu_s)

            update_progress(80, "Calcolo radianza up-welling...")
            # Calcolo della radianza up-welling (combinando superficie ed atmosfera)
            radianza_upwelling = surface_emissivity * blackbodyn(surface_temperature, w) * tau[:, 0]
            for l in range(n_layers - 1, -1, -1):
                dtau = tau[:, l+1] - tau[:, l]  # variazione discreta della trasmittanza
                radianza_upwelling += blackbodyn(tm[l], w) * dtau

            update_progress(85, "Calcolo radianza down-welling...")
            # Calcolo della radianza down-welling atmosferica
            radianza_downwelling = np.zeros_like(w)
            for l in range(n_layers):
                dtau_inv = 1/tau[:, l+1] - 1/tau[:, l]  # variazione discreta
                radianza_downwelling += blackbodyn(tm[l], w) * dtau_inv
            radianza_downwelling *= (surface_emissivity - 1) * tau[:, 0]**2

            update_progress(90, "Calcolo contributo solare...")
            # Calcolo del contributo solare diretto
            T_sole = 5800
            dim_ang_sole = 6.8e-5
            E_sigma = blackbodyn(T_sole, w) * dim_ang_sole

            radianza_solare = (1 - surface_emissivity)/np.pi * tau[:, 0] * \
                               tau_mu[:, 0] * mu_s * E_sigma

            update_progress(95, "Calcolo radianza totale...")
            # Calcolo della radianza totale includendo tutte le componenti
            radianza_totale = (radianza_upwelling + 
                               radianza_downwelling + 
                               radianza_solare)

            # Chiude la finestra di progresso
            progress_window.destroy()

            # Display radiance plots in a new pop-up window with themed styling
            radiance_window = tk.Toplevel(self)
            radiance_window.title("Grafici delle Radianze")
            radiance_window.geometry("800x600")
            radiance_window.configure(bg=self.get_current_color('bg'))

            # Crea un canvas scrollabile
            canvas = tk.Canvas(radiance_window, bg=self.get_current_color('bg'))
            scrollbar_y = tk.Scrollbar(radiance_window, orient=tk.VERTICAL, command=canvas.yview, 
                                       bg=self.get_current_color('button_bg'), 
                                       troughcolor=self.get_current_color('frame_bg'))
            
            # Configura il canvas
            canvas.configure(yscrollcommand=scrollbar_y.set)
            
            # Posiziona lo scrollbar
            scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)
            canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            
            # Crea un frame interno al canvas per contenere i grafici
            inner_frame = tk.Frame(canvas, bg=self.get_current_color('bg'))
            canvas.create_window((0, 0), window=inner_frame, anchor='nw', width=radiance_window.winfo_width())
            
            # Creazione delle figure con dimensione iniziale fissa
            fig1 = Figure(figsize=(8, 4), dpi=100, facecolor=self.get_current_color('bg'))
            ax1 = fig1.add_subplot(111)
            ax1.plot(w, radianza_upwelling, 'b', linewidth=0.8, label='Up-welling')
            ax1.plot(w, radianza_downwelling, 'r', linewidth=0.8, label='Down-welling')
            ax1.plot(w, radianza_solare, 'g', linewidth=0.8, label='Contributo Solare Diretto')
            ax1.legend()
            ax1.set_xlabel("Numero d'onda (cm$^{-1}$)")
            ax1.set_ylabel("Radianza")
            ax1.set_title("Componenti Radiative Dirette")
            ax1.grid(True)
            fig1.tight_layout(pad=0)
            
            fig2 = Figure(figsize=(8, 4), dpi=100, facecolor=self.get_current_color('bg'))
            ax2 = fig2.add_subplot(111)
            ax2.plot(w, radianza_solare, 'c', linewidth=0.8, label='Solare Riflessa')
            ax2.plot(w, radianza_downwelling, 'm', linewidth=0.8, label='Atmosferica Riflessa')
            ax2.legend()
            ax2.set_xlabel("Numero d'onda (cm$^{-1}$)")
            ax2.set_ylabel("Radianza")
            ax2.set_title("Componenti Riflesse")
            ax2.grid(True)
            fig2.tight_layout(pad=0)
            
            fig3 = Figure(figsize=(8, 4), dpi=100, facecolor=self.get_current_color('bg'))
            ax3 = fig3.add_subplot(111)
            ax3.plot(w, radianza_totale, 'k', linewidth=0.8)
            ax3.set_xlabel("Numero d'onda (cm$^{-1}$)")
            ax3.set_ylabel("Radianza")
            ax3.set_title("Radianza Totale (Tutte le Componenti)")
            ax3.grid(True)
            fig3.tight_layout(pad=0)
            
            # Personalizza i colori degli assi per adattarsi al tema
            for ax in [ax1, ax2, ax3]:
                ax.set_facecolor(self.get_current_color('frame_bg'))
                ax.spines['bottom'].set_color(self.get_current_color('text_primary'))
                ax.spines['top'].set_color(self.get_current_color('text_primary'))
                ax.spines['right'].set_color(self.get_current_color('text_primary'))
                ax.spines['left'].set_color(self.get_current_color('text_primary'))
                ax.tick_params(colors=self.get_current_color('text_primary'))
                ax.xaxis.label.set_color(self.get_current_color('text_primary'))
                ax.yaxis.label.set_color(self.get_current_color('text_primary'))
                ax.title.set_color(self.get_current_color('text_primary'))
            
            # Crea i canvas per ogni figura
            canvas_plot1 = FigureCanvasTkAgg(fig1, master=inner_frame)
            canvas_widget1 = canvas_plot1.get_tk_widget()
            canvas_widget1.configure(bg=self.get_current_color('bg'))
            canvas_widget1.pack(fill=tk.BOTH, expand=True)
            
            canvas_plot2 = FigureCanvasTkAgg(fig2, master=inner_frame)
            canvas_widget2 = canvas_plot2.get_tk_widget()
            canvas_widget2.configure(bg=self.get_current_color('bg'))
            canvas_widget2.pack(fill=tk.BOTH, expand=True)
            
            canvas_plot3 = FigureCanvasTkAgg(fig3, master=inner_frame)
            canvas_widget3 = canvas_plot3.get_tk_widget()
            canvas_widget3.configure(bg=self.get_current_color('bg'))
            canvas_widget3.pack(fill=tk.BOTH, expand=True)
            
            # Aggiungi le barre degli strumenti
            toolbar1 = NavigationToolbar2Tk(canvas_plot1, inner_frame)
            toolbar1.update()
            toolbar1.pack(fill=tk.X)
            
            toolbar2 = NavigationToolbar2Tk(canvas_plot2, inner_frame)
            toolbar2.update()
            toolbar2.pack(fill=tk.X)
            
            toolbar3 = NavigationToolbar2Tk(canvas_plot3, inner_frame)
            toolbar3.update()
            toolbar3.pack(fill=tk.X)
            
            # Funzione per aggiornare dinamicamente le dimensioni dei grafici
            def on_resize(event):
                # Aggiorna la larghezza del frame interno
                canvas.itemconfig("all", width=event.width)
                
                # Aggiorna le dimensioni delle figure
                new_width = event.width / 100  # Converti da pixel a pollici
                new_height = event.height / 100
                
                for fig, canvas_plot in [(fig1, canvas_plot1), (fig2, canvas_plot2), (fig3, canvas_plot3)]:
                    fig.set_size_inches(new_width, new_height, forward=True)
                    canvas_plot.draw()
            
            # Associa il binding di resize al canvas principale
            canvas.bind("<Configure>", on_resize)
            
            # Aggiorna la regione scrollabile
            inner_frame.update_idletasks()
            canvas.config(scrollregion=canvas.bbox("all"))
            
            # Disegna i canvas
            canvas_plot1.draw()
            canvas_plot2.draw()
            canvas_plot3.draw()

        except Exception as e:
            if 'progress_window' in locals():
                progress_window.destroy()
            messagebox.showerror("Errore", f"Errore nel calcolo della radianza: {e}")

# =============================================================================
# PROGRAMMA PRINCIPALE
# =============================================================================

if __name__ == "__main__":
    app = RadianceGUI()
    app.mainloop()
