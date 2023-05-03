/* TO COMPILE:
	g++ -fopenmp -c main.cpp -I../lib/SFML-2.5.1/include
	g++ main.o -fopenmp -o app -L../lib/SFML-2.5.1/lib -lsfml-graphics -lsfml-window -lsfml-system
   TO RUN:
	export LD_LIBRARY_PATH=../lib/SFML-2.5.1/lib && ./app
*/

#include <SFML/Graphics.hpp>
#include "FPS.cpp"
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <omp.h>
#include <iostream>

struct Body
{
    sf::Vector3f pos;
    sf::Vector3f vel; // 2 m/s along x-axis
    sf::Vector3f acc; // no acceleration at first
    float mass = 1.0f; // 1kg
    float drag = 1.0f; // rho*C*Area � simplified drag for this example
    sf::CircleShape shape;

    Body() {
        this->pos = sf::Vector3f(0.0f, 0.0f, 0.0f);
        this->vel = sf::Vector3f(50.0f, 0.0f, 0.0f);
        this->acc = sf::Vector3f(0.0f, 0.0f, 0.0f);
        this->shape = sf::CircleShape(10.f);
        this->shape.setOrigin(this->shape.getRadius(), this->shape.getRadius());
        this->shape.setPosition(this->pos.x, this->pos.y);
    }

    /**
     * Update pos and vel using "Velocity Verlet" integration
     * @param dt DeltaTime / time step [eg: 0.01]
     */
    void update(float dt)
    {
        sf::Vector3f new_pos = this->pos + this->vel * dt + this->acc * (dt * dt * 0.5f);
        sf::Vector3f new_acc = apply_forces(); // only needed if acceleration is not constant
        //sf::Vector3f new_acc = sf::Vector3f(0.0f, 9.81f, 0.0f);
        sf::Vector3f new_vel = this->vel + (this->acc + new_acc) * (dt * 0.5f);
        this->pos = new_pos;
        this->vel = new_vel;
        this->acc = new_acc;
        this->shape.setPosition(this->pos.x, this->pos.y);
    }

    void setPosition(float x, float y) 
    {
        this->pos = sf::Vector3f(x, y, 0.0f);
        this->shape.setPosition(this->pos.x, this->pos.y);
    }

    void updatePositionByShape()
    {
        this->pos = sf::Vector3f(this->shape.getPosition().x, this->shape.getPosition().y, 0.0f);
    }

    void updatePositionByBody()
    {
        this->shape.setPosition(this->pos.x, this->pos.y);
    }

    sf::Vector3f apply_forces() const
    {
        sf::Vector3f grav_acc = sf::Vector3f( 0.0f, 40.0f, 0.0f ); // 9.81 m/s� down in the y-axis
        // sf::Vector3f newVel = sf::Vector3f(pow(vel.x,2.0f), pow(vel.y, 2.0f), pow(vel.z, 2.0f));
        sf::Vector3f drag_force = 0.5f * drag * vel; // D = 0.5 * (rho * C * Area * vel^2)
        sf::Vector3f drag_acc = drag_force / mass; // a = F/m
        return grav_acc- drag_acc;
    }
};

static const size_t screenX = 800;
static const size_t screenY = 800;
std::vector<Body*> particles;
void findCollisions();
void updatePhysics(float dt);
void updatePhysicsSubtick(float dt, int subTicks);

int main()
{
    omp_set_num_threads(8);
    sf::RenderWindow window(sf::VideoMode(800, 800), "SFML works!");
    window.setFramerateLimit(60);
    const float dt = 1.0f / static_cast<float>(60);
    particles.push_back(new Body());

    FPS fps;
    unsigned int counter = 0;
    int colorCounter = 0;

    sf::Font font; //set font
    if(!font.loadFromFile("/usr/share/fonts/truetype/ubuntu/UbuntuMono-BI.ttf"))
    {
        std::cout << "Error loading font!" << std::endl;
    }

    sf::Text text;
    text.setFont(font);
    text.setCharacterSize(20);
    text.setFillColor(sf::Color::White);
    text.setPosition(150.f, 0.f);

    sf::Clock clock; //start clock

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed){
				window.close();
			}
                
        }
        if (counter % 15 == 0 && particles.size() < 100)
        {
            Body* newBody = new Body();
            if (colorCounter == 0) newBody->shape.setFillColor(sf::Color(255, 0, 0, 255));
            if (colorCounter == 1) newBody->shape.setFillColor(sf::Color(0, 255, 0, 255));
            if (colorCounter == 2) newBody->shape.setFillColor(sf::Color(0, 0, 255, 255));
            if (colorCounter < 2) colorCounter++;
            else colorCounter = 0;
            particles.push_back(newBody);
        }

        updatePhysicsSubtick(dt, 8);

        fps.update();
        std::ostringstream ss;
        ss << fps.getFPS();

        window.setTitle(ss.str());

        window.clear();
        for (Body* particle : particles)
        {
            window.draw(particle->shape);
        }

        //display time and num particles
        sf::Time elapsed = clock.getElapsedTime();
        text.setString("Elapsed time: " + std::to_string(elapsed.asSeconds()) + "    Number of Particles: " + std::to_string(particles.size()));
        window.draw(text);

        window.display();
        if (counter < 60) counter++;
        else counter = 0;
    }

    sf::Time elapsed = clock.getElapsedTime();
    std::cout << "Elapsed time: " << elapsed.asSeconds() << std::endl;

    particles.clear();
    return 0;
}

float distance(int x1, int y1, int x2, int y2)
{
    // Calculating distance
    return (float)sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}

bool collide(Body* particle1, Body* particle2)
{
    if((float)distance(particle1->pos.x, particle1->pos.y, particle2->pos.x, particle2->pos.y) < (particle1->shape.getRadius() + particle2->shape.getRadius()))
        return true;
    else
        return false;
}

void solveCollision(Body* particle1, Body* particle2)
{
    constexpr float response_coef = 1.0f;
    const sf::Vector2f o2_o1 = particle1->shape.getPosition() - particle2->shape.getPosition();
    const float dist2 = o2_o1.x * o2_o1.x + o2_o1.y * o2_o1.y;
    
    const float dist = sqrt(dist2);
    const float delta = response_coef * 0.5f * ((particle1->shape.getRadius()+particle2->shape.getRadius()) - dist);
    //const float delta = dist / 4.0f;
    const sf::Vector2f col_vec = (o2_o1 / dist) * delta;
    sf::Vector2f newP1Pos = particle1->shape.getPosition() + col_vec;
    particle1->shape.setPosition(newP1Pos);
    particle1->updatePositionByShape();
    sf::Vector2f newP2Pos = particle2->shape.getPosition() - col_vec;
    particle2->shape.setPosition(newP2Pos);
    particle2->updatePositionByShape();
}

void findCollisions() 
{
	int i;
	#pragma omp parallel for
    for (i =0; i < particles.size(); i++)
    {
		auto particle1 = particles.at(i);
        for (auto particle2 : particles)
        {
            if (particle1 != particle2)
            {
                if (collide(particle1, particle2))
                {
                    solveCollision(particle1, particle2);
                }
            }
        }
    }
}

void updatePhysics(float dt)
{
    const float margin = 2.0f;
    int i;

	// fillCells();
	// find_collisions_grid();
    findCollisions();

    #pragma omp parallel for
    for (i =0; i < particles.size(); i++)
    {
		auto particle = particles.at(i);
		
        //findCollisions();
        
        if (particle->pos.x > 800 - margin - particle->shape.getRadius()) {
            particle->pos.x = 800 - margin - particle->shape.getRadius();
            particle->updatePositionByBody();
			particle->vel = sf::Vector3f(-particle->vel.x,particle->vel.y,particle->vel.z);
        }
        else if (particle->pos.x < margin + particle->shape.getRadius()) {
            particle->pos.x = margin + particle->shape.getRadius();
            particle->updatePositionByBody();
			particle->vel = sf::Vector3f(-particle->vel.x,particle->vel.y,particle->vel.z);
        }
        if (particle->pos.y > 800 - margin - particle->shape.getRadius()) {
            particle->pos.y = 800 - margin - particle->shape.getRadius();
            particle->updatePositionByBody();
			particle->vel = sf::Vector3f(particle->vel.x,-particle->vel.y,particle->vel.z);
        }
        else if (particle->pos.y < margin + particle->shape.getRadius()) {
            particle->pos.y = margin + particle->shape.getRadius();
            particle->updatePositionByBody();
			particle->vel = sf::Vector3f(particle->vel.x,-particle->vel.y,particle->vel.z);
        }
		particle->update(dt);
    }
}

void updatePhysicsSubtick(float dt, int subTicks)
{
    const float sub_dt = dt / (float)subTicks;
    for (int i{ subTicks }; i--;)
    {
        updatePhysics(sub_dt);
    }
}