# N-Body Simulation in C++ with OpenGL

Diese Anwendung ist eine interaktive 3D-N-Body-Simulation, geschrieben in C++ und visualisiert mit OpenGL. Sie simuliert gravitative Wechselwirkungen zwischen mehreren Partikeln in einem dreidimensionalen Raum.

## 🧠 Features

- 3D-Visualisierung von Partikelsystemen unter Gravitation
- Konfigurierbare Parameter wie:
  - Anzahl der Partikel
  - Simulationsraumgröße
  - Zeitauflösung
  - Zentrale Masse (optional)
- Verschiedene Integrationsmethoden:
  - Euler
  - Verlet
  - Runge-Kutta (RK4)
- Spurenanzeige für Bewegungspfad der Partikel
- Unterstützung für macOS und Linux

## 🖥️ Screenshots

![Simulation Vorschau](Screenshot.png)

## 🔧 Voraussetzungen

- C++17 oder neuer
- OpenGL
- GLUT (z. B. freeglut oder GLUT.framework auf macOS)

### macOS

```bash
brew install freeglut
```

### Linux (Debian/Ubuntu)

```bash
sudo apt install freeglut3-dev
```

## 🧪 Kompilierung

### macOS:

```bash
g++ main.cpp -o NBodySimulation -framework OpenGL -framework GLUT
```

### Linux:

```bash
g++ main.cpp -o NBodySimulation -lGL -lGLU -lglut
```

## 🚀 Ausführung

```bash
./NBodySimulation
```

## ⚙️ Konfigurierbare Parameter

Die Parameter sind in der Struktur `SimulationConfig` im Quellcode definiert:

```cpp
struct SimulationConfig {
    double timeStep = 0.01;
    double softening = 0.1;
    int particleCount = 100;
    double spaceSize = 10.0;
    double maxInitialVelocity = 0.5;
    double centralMass = 1000.0;
    bool useCentralMass = true;
    bool showTrails = true;
    int trailLength = 50;
    std::string integrationMethod = "verlet"; // Optionen: euler, verlet, rk4
};
```

## 📁 Projektstruktur

```
main.cpp        // Hauptprogramm mit OpenGL-Rendering und Simulation
```

## 🧑‍💻 Entwickler

**Bengin Sternas**  
Erstellt am 15.07.2025

## 📜 Lizenz

Dieses Projekt steht unter der MIT-Lizenz. Weitere Informationen findest du in der Datei `LICENSE`
