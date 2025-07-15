//
//  main.cpp
//  NBodySimulation
//
//  Created by Bengin Sternas on 15.07.25.
//

/**
 * N-Body Gravity Simulation
 *
 * Eine komplexe Simulation der gravitativen Wechselwirkungen zwischen mehreren Körpern
 * unter Verwendung der Newton'schen Gravitationsgesetze.
 *
 * Features:
 * - Objektorientierte Architektur
 * - OpenGL-basierte Visualisierung
 * - Konfigurierbare Parameter
 * - Verschiedene Integrationsmethoden (Euler, Verlet, RK4)
 * - Energie- und Impulserhaltung-Monitoring
 *
 * Kompilierung in Xcode:
 * 1. Neues Command Line Tool Projekt erstellen (C++)
 * 2. OpenGL und GLUT Frameworks hinzufügen
 * 3. Build Settings -> Other Linker Flags: -framework OpenGL -framework GLUT
 * 4. Für Warnungen: GL_SILENCE_DEPRECATION definieren (siehe unten)
 */

// Unterdrückt OpenGL Deprecation Warnungen
#define GL_SILENCE_DEPRECATION

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <sstream>
#include <memory>
#include <string>
#include <algorithm>

// Mathematische Konstanten
const double PI = 3.14159265358979323846;
const double G = 6.67430e-11; // Gravitationskonstante (skaliert für Simulation)
const double SIMULATION_G = 1.0; // Angepasste Gravitationskonstante für bessere Visualisierung

// Simulationsparameter
struct SimulationConfig {
    double timeStep = 0.01;
    double softening = 0.1; // Vermeidet Singularitäten bei kleinen Abständen
    int particleCount = 100;
    double spaceSize = 10.0;
    double maxInitialVelocity = 0.5;
    double centralMass = 1000.0; // Masse des zentralen Körpers
    bool useCentralMass = true;
    bool showTrails = true;
    int trailLength = 50;
    std::string integrationMethod = "verlet"; // euler, verlet, rk4
};

// 3D Vektor Klasse
class Vector3D {
public:
    double x, y, z;
    
    Vector3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    
    Vector3D operator+(const Vector3D& v) const {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }
    
    Vector3D operator-(const Vector3D& v) const {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }
    
    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }
    
    Vector3D operator/(double scalar) const {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }
    
    Vector3D& operator+=(const Vector3D& v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    
    double magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }
    
    double magnitudeSquared() const {
        return x * x + y * y + z * z;
    }
    
    Vector3D normalized() const {
        double mag = magnitude();
        if (mag == 0) return Vector3D(0, 0, 0);
        return *this / mag;
    }
    
    double dot(const Vector3D& v) const {
        return x * v.x + y * v.y + z * v.z;
    }
    
    Vector3D cross(const Vector3D& v) const {
        return Vector3D(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }
};

// Körper-Klasse (repräsentiert einen Himmelskörper)
class Body {
private:
    std::vector<Vector3D> trail;
    int maxTrailLength;
    
public:
    Vector3D position;
    Vector3D velocity;
    Vector3D acceleration;
    Vector3D previousAcceleration; // Für Verlet-Integration
    double mass;
    double radius;
    Vector3D color;
    bool isFixed; // Für zentrale Masse
    
    Body(const Vector3D& pos, const Vector3D& vel, double m, double r = 0.1)
        : position(pos), velocity(vel), mass(m), radius(r),
          acceleration(0, 0, 0), previousAcceleration(0, 0, 0),
          isFixed(false), maxTrailLength(50) {
        // Zufällige Farbe generieren
        color = Vector3D(
            0.3 + (rand() % 70) / 100.0,
            0.3 + (rand() % 70) / 100.0,
            0.3 + (rand() % 70) / 100.0
        );
    }
    
    void updateTrail() {
        trail.push_back(position);
        if (static_cast<int>(trail.size()) > maxTrailLength) {
            trail.erase(trail.begin());
        }
    }
    
    void setTrailLength(int length) {
        maxTrailLength = length;
        while (static_cast<int>(trail.size()) > maxTrailLength && !trail.empty()) {
            trail.erase(trail.begin());
        }
    }
    
    const std::vector<Vector3D>& getTrail() const {
        return trail;
    }
    
    void clearTrail() {
        trail.clear();
    }
};

// Integrator-Klasse für verschiedene numerische Integrationsmethoden
class Integrator {
public:
    virtual ~Integrator() = default;
    virtual void integrate(std::vector<Body>& bodies, double dt, double softening) = 0;
    
protected:
    Vector3D calculateAcceleration(const Body& body, const std::vector<Body>& bodies,
                                   size_t bodyIndex, double softening) {
        Vector3D acc(0, 0, 0);
        
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (i == bodyIndex) continue;
            
            Vector3D r = bodies[i].position - body.position;
            double distSq = r.magnitudeSquared() + softening * softening;
            double dist = std::sqrt(distSq);
            double force = SIMULATION_G * bodies[i].mass / (distSq * dist);
            
            acc += r * force;
        }
        
        return acc;
    }
};

// Euler-Integration
class EulerIntegrator : public Integrator {
public:
    void integrate(std::vector<Body>& bodies, double dt, double softening) override {
        // Beschleunigungen berechnen
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (!bodies[i].isFixed) {
                bodies[i].acceleration = calculateAcceleration(bodies[i], bodies, i, softening);
            }
        }
        
        // Positionen und Geschwindigkeiten aktualisieren
        for (auto& body : bodies) {
            if (!body.isFixed) {
                body.velocity += body.acceleration * dt;
                body.position += body.velocity * dt;
            }
        }
    }
};

// Verlet-Integration (bessere Energieerhaltung)
class VerletIntegrator : public Integrator {
public:
    void integrate(std::vector<Body>& bodies, double dt, double softening) override {
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (bodies[i].isFixed) continue;
            
            // Neue Beschleunigung berechnen
            Vector3D newAcceleration = calculateAcceleration(bodies[i], bodies, i, softening);
            
            // Verlet-Integration
            bodies[i].position += bodies[i].velocity * dt +
                                  bodies[i].acceleration * (0.5 * dt * dt);
            bodies[i].velocity += (bodies[i].acceleration + newAcceleration) * (0.5 * dt);
            
            bodies[i].previousAcceleration = bodies[i].acceleration;
            bodies[i].acceleration = newAcceleration;
        }
    }
};

// Runge-Kutta 4. Ordnung (RK4) - höchste Genauigkeit
class RK4Integrator : public Integrator {
private:
    struct State {
        Vector3D position;
        Vector3D velocity;
    };
    
    struct Derivative {
        Vector3D velocity;
        Vector3D acceleration;
    };
    
    Derivative evaluate(const Body& body, const std::vector<Body>& bodies,
                       size_t bodyIndex, double softening, double dt,
                       const Derivative& d) {
        State state;
        state.position = body.position + d.velocity * dt;
        state.velocity = body.velocity + d.acceleration * dt;
        
        Body tempBody = body;
        tempBody.position = state.position;
        tempBody.velocity = state.velocity;
        
        Derivative output;
        output.velocity = state.velocity;
        output.acceleration = calculateAcceleration(tempBody, bodies, bodyIndex, softening);
        
        return output;
    }
    
public:
    void integrate(std::vector<Body>& bodies, double dt, double softening) override {
        std::vector<State> newStates(bodies.size());
        
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (bodies[i].isFixed) {
                newStates[i].position = bodies[i].position;
                newStates[i].velocity = bodies[i].velocity;
                continue;
            }
            
            Derivative k1, k2, k3, k4;
            
            k1.velocity = bodies[i].velocity;
            k1.acceleration = calculateAcceleration(bodies[i], bodies, i, softening);
            
            k2 = evaluate(bodies[i], bodies, i, softening, dt * 0.5, k1);
            k3 = evaluate(bodies[i], bodies, i, softening, dt * 0.5, k2);
            k4 = evaluate(bodies[i], bodies, i, softening, dt, k3);
            
            Vector3D dxdt = (k1.velocity + k2.velocity * 2.0 + k3.velocity * 2.0 + k4.velocity) / 6.0;
            Vector3D dvdt = (k1.acceleration + k2.acceleration * 2.0 + k3.acceleration * 2.0 + k4.acceleration) / 6.0;
            
            newStates[i].position = bodies[i].position + dxdt * dt;
            newStates[i].velocity = bodies[i].velocity + dvdt * dt;
        }
        
        // Neue Zustände übernehmen
        for (size_t i = 0; i < bodies.size(); ++i) {
            bodies[i].position = newStates[i].position;
            bodies[i].velocity = newStates[i].velocity;
        }
    }
};

// Haupt-Simulationsklasse
class NBodySimulation {
private:
    std::vector<Body> bodies;
    SimulationConfig config;
    std::unique_ptr<Integrator> integrator;
    double totalTime;
    int frameCount;
    
    // Kamera-Parameter
    float cameraDistance;
    float cameraRotationX;
    float cameraRotationY;
    
    // Statistiken
    double totalEnergy;
    Vector3D totalMomentum;
    
public:
    NBodySimulation(const SimulationConfig& cfg)
        : config(cfg), totalTime(0), frameCount(0),
          cameraDistance(30.0f), cameraRotationX(30.0f), cameraRotationY(0.0f) {
        
        initializeIntegrator();
        initializeBodies();
    }
    
    void initializeIntegrator() {
        if (config.integrationMethod == "euler") {
            integrator = std::make_unique<EulerIntegrator>();
        } else if (config.integrationMethod == "rk4") {
            integrator = std::make_unique<RK4Integrator>();
        } else {
            integrator = std::make_unique<VerletIntegrator>();
        }
    }
    
    void initializeBodies() {
        bodies.clear();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1.0, 1.0);
        std::uniform_real_distribution<> massDis(0.5, 2.0);
        
        // Zentraler massereicher Körper (optional)
        if (config.useCentralMass) {
            Body centralBody(Vector3D(0, 0, 0), Vector3D(0, 0, 0),
                           config.centralMass, 0.5);
            centralBody.isFixed = true;
            centralBody.color = Vector3D(1.0, 0.8, 0.0); // Gelb für Sonne
            bodies.push_back(centralBody);
        }
        
        // Andere Körper
        for (int i = 0; i < config.particleCount; ++i) {
            // Zufällige Position in einer Scheibe
            double angle = 2 * PI * dis(gen);
            double radius = config.spaceSize * (0.3 + 0.7 * std::abs(dis(gen)));
            double height = config.spaceSize * 0.1 * dis(gen); // Flache Scheibe
            
            Vector3D pos(radius * std::cos(angle), height, radius * std::sin(angle));
            
            // Orbital-Geschwindigkeit für stabilen Orbit
            double orbitalSpeed = 0;
            if (config.useCentralMass) {
                orbitalSpeed = std::sqrt(SIMULATION_G * config.centralMass / radius);
            }
            
            // Geschwindigkeit tangential zur Umlaufbahn
            Vector3D vel(-orbitalSpeed * std::sin(angle), 0, orbitalSpeed * std::cos(angle));
            
            // Kleine zufällige Störung hinzufügen
            vel.x += config.maxInitialVelocity * 0.1 * dis(gen);
            vel.y += config.maxInitialVelocity * 0.1 * dis(gen);
            vel.z += config.maxInitialVelocity * 0.1 * dis(gen);
            
            double mass = massDis(gen);
            Body body(pos, vel, mass, 0.05 + 0.05 * mass);
            bodies.push_back(body);
        }
        
        updateTrailLength();
    }
    
    void update() {
        integrator->integrate(bodies, config.timeStep, config.softening);
        
        totalTime += config.timeStep;
        frameCount++;
        
        // Trails aktualisieren
        if (config.showTrails && frameCount % 2 == 0) {
            for (auto& body : bodies) {
                body.updateTrail();
            }
        }
        
        // Statistiken alle 30 Frames berechnen
        if (frameCount % 30 == 0) {
            calculateStatistics();
        }
    }
    
    void calculateStatistics() {
        totalEnergy = 0;
        totalMomentum = Vector3D(0, 0, 0);
        
        // Kinetische Energie und Impuls
        for (const auto& body : bodies) {
            double kineticEnergy = 0.5 * body.mass * body.velocity.magnitudeSquared();
            totalEnergy += kineticEnergy;
            totalMomentum += body.velocity * body.mass;
        }
        
        // Potentielle Energie
        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                Vector3D r = bodies[j].position - bodies[i].position;
                double dist = r.magnitude();
                if (dist > 0) {
                    double potentialEnergy = -SIMULATION_G * bodies[i].mass * bodies[j].mass / dist;
                    totalEnergy += potentialEnergy;
                }
            }
        }
    }
    
    void render() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        // Kamera-Position
        glTranslatef(0, 0, -cameraDistance);
        glRotatef(cameraRotationX, 1, 0, 0);
        glRotatef(cameraRotationY, 0, 1, 0);
        
        // Koordinatenachsen
        drawAxes();
        
        // Körper und Trails rendern
        for (const auto& body : bodies) {
            // Trail
            if (config.showTrails) {
                drawTrail(body);
            }
            
            // Körper
            glColor3f(static_cast<float>(body.color.x),
                     static_cast<float>(body.color.y),
                     static_cast<float>(body.color.z));
            glPushMatrix();
            glTranslatef(static_cast<float>(body.position.x),
                        static_cast<float>(body.position.y),
                        static_cast<float>(body.position.z));
            glutSolidSphere(body.radius, 20, 20);
            glPopMatrix();
        }
        
        // Statistiken anzeigen
        displayStatistics();
        
        glutSwapBuffers();
    }
    
    void drawTrail(const Body& body) {
        const auto& trail = body.getTrail();
        if (trail.size() < 2) return;
        
        glBegin(GL_LINE_STRIP);
        for (size_t i = 0; i < trail.size(); ++i) {
            float alpha = static_cast<float>(i) / static_cast<float>(trail.size());
            glColor4f(static_cast<float>(body.color.x),
                     static_cast<float>(body.color.y),
                     static_cast<float>(body.color.z),
                     alpha * 0.5f);
            glVertex3f(static_cast<float>(trail[i].x),
                      static_cast<float>(trail[i].y),
                      static_cast<float>(trail[i].z));
        }
        glEnd();
    }
    
    void drawAxes() {
        glLineWidth(2.0f);
        
        // X-Achse (rot)
        glColor3f(1, 0, 0);
        glBegin(GL_LINES);
        glVertex3f(0, 0, 0);
        glVertex3f(5, 0, 0);
        glEnd();
        
        // Y-Achse (grün)
        glColor3f(0, 1, 0);
        glBegin(GL_LINES);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 5, 0);
        glEnd();
        
        // Z-Achse (blau)
        glColor3f(0, 0, 1);
        glBegin(GL_LINES);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 0, 5);
        glEnd();
        
        glLineWidth(1.0f);
    }
    
    void displayStatistics() {
        // 2D-Modus für Text
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, 800, 600, 0, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        
        glColor3f(1, 1, 1);
        
        // Statistiken formatieren
        std::vector<std::string> stats = {
            "N-Body Gravity Simulation",
            "Bodies: " + std::to_string(bodies.size()),
            "Time: " + std::to_string(static_cast<int>(totalTime)),
            "Energy: " + std::to_string(static_cast<int>(totalEnergy)),
            "Integration: " + config.integrationMethod,
            "",
            "Controls:",
            "R - Reset",
            "Space - Pause",
            "T - Toggle trails",
            "+/- Change speed",
            "1/2/3 - Integration method"
        };
        
        int y = 20;
        for (const auto& stat : stats) {
            glRasterPos2i(10, y);
            for (char c : stat) {
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, c);
            }
            y += 15;
        }
        
        // 3D-Modus wiederherstellen
        glPopMatrix();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
    }
    
    // Getter und Setter
    void setCameraRotation(float x, float y) {
        cameraRotationX = x;
        cameraRotationY = y;
    }
    
    void setCameraDistance(float dist) {
        cameraDistance = std::max(5.0f, std::min(100.0f, dist));
    }
    
    void changeTimeStep(double factor) {
        config.timeStep *= factor;
        config.timeStep = std::max(0.0001, std::min(0.1, config.timeStep));
    }
    
    void toggleTrails() {
        config.showTrails = !config.showTrails;
        if (!config.showTrails) {
            for (auto& body : bodies) {
                body.clearTrail();
            }
        }
    }
    
    void setIntegrationMethod(const std::string& method) {
        config.integrationMethod = method;
        initializeIntegrator();
    }
    
    void updateTrailLength() {
        for (auto& body : bodies) {
            body.setTrailLength(config.trailLength);
        }
    }
    
    void reset() {
        initializeBodies();
        totalTime = 0;
        frameCount = 0;
    }
    
    float getCameraRotationX() const { return cameraRotationX; }
    float getCameraRotationY() const { return cameraRotationY; }
    float getCameraDistance() const { return cameraDistance; }
};

// Globale Variablen für GLUT-Callbacks
NBodySimulation* g_simulation = nullptr;
bool g_paused = false;
int g_mouseX = 0, g_mouseY = 0;
bool g_mouseLeftDown = false;
bool g_mouseRightDown = false;

// GLUT-Callbacks
void display() {
    if (g_simulation) {
        g_simulation->render();
    }
}

void idle() {
    if (g_simulation && !g_paused) {
        g_simulation->update();
        glutPostRedisplay();
    }
}

void reshape(int width, int height) {
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, static_cast<double>(width) / height, 0.1, 1000.0);
    glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y) {
    if (!g_simulation) return;
    
    switch (key) {
        case 'r':
        case 'R':
            g_simulation->reset();
            break;
        case ' ':
            g_paused = !g_paused;
            break;
        case 't':
        case 'T':
            g_simulation->toggleTrails();
            break;
        case '+':
        case '=':
            g_simulation->changeTimeStep(1.2);
            break;
        case '-':
        case '_':
            g_simulation->changeTimeStep(0.8);
            break;
        case '1':
            g_simulation->setIntegrationMethod("euler");
            break;
        case '2':
            g_simulation->setIntegrationMethod("verlet");
            break;
        case '3':
            g_simulation->setIntegrationMethod("rk4");
            break;
        case 'q':
        case 'Q':
        case 27: // ESC
            delete g_simulation;
            exit(0);
            break;
    }
    
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
    g_mouseX = x;
    g_mouseY = y;
    
    if (button == GLUT_LEFT_BUTTON) {
        g_mouseLeftDown = (state == GLUT_DOWN);
    } else if (button == GLUT_RIGHT_BUTTON) {
        g_mouseRightDown = (state == GLUT_DOWN);
    }
}

void motion(int x, int y) {
    if (!g_simulation) return;
    
    int dx = x - g_mouseX;
    int dy = y - g_mouseY;
    
    if (g_mouseLeftDown) {
        float rotY = g_simulation->getCameraRotationY() + dx * 0.5f;
        float rotX = g_simulation->getCameraRotationX() + dy * 0.5f;
        g_simulation->setCameraRotation(rotX, rotY);
    } else if (g_mouseRightDown) {
        float dist = g_simulation->getCameraDistance() + dy * 0.1f;
        g_simulation->setCameraDistance(dist);
    }
    
    g_mouseX = x;
    g_mouseY = y;
    glutPostRedisplay();
}

// Hauptfunktion
int main(int argc, char** argv) {
    // GLUT initialisieren
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("N-Body Gravity Simulation");
    
    // OpenGL-Einstellungen
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    
    // Simulation erstellen
    SimulationConfig config;
    config.particleCount = 100;
    config.useCentralMass = true;
    config.showTrails = true;
    config.integrationMethod = "verlet";
    
    g_simulation = new NBodySimulation(config);
    
    // GLUT-Callbacks registrieren
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    
    // Hauptschleife starten
    std::cout << "N-Body Gravity Simulation gestartet!" << std::endl;
    std::cout << "Steuerung:" << std::endl;
    std::cout << "  R - Reset" << std::endl;
    std::cout << "  Leertaste - Pause" << std::endl;
    std::cout << "  T - Trails ein/aus" << std::endl;
    std::cout << "  +/- - Geschwindigkeit ändern" << std::endl;
    std::cout << "  1/2/3 - Integrationsmethode wechseln" << std::endl;
    std::cout << "  Maus - Kamera drehen/zoomen" << std::endl;
    std::cout << "  ESC/Q - Beenden" << std::endl;
    
    glutMainLoop();
    
    return 0;
}
