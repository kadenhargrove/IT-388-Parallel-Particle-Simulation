/* TO COMPILE:
	g++ -fopenmp -g -c main.cpp -I../lib/SFML-2.5.1/include
	g++ main.o -fopenmp -o app -L../lib/SFML-2.5.1/lib -lsfml-graphics -lsfml-window -lsfml-system
   TO RUN:
	export LD_LIBRARY_PATH=../lib/SFML-2.5.1/lib && ./app nThreads nParticles cellSize drawWindowFlag
    cellSize must be divisible by screen x and y size (800) and >= 2x particle radius (5.0), drawWindowFlag can be 0 or 1

   DEBUG: g++ -fopenmp -g -Wall -Wextra -pedantic -c main.cpp -I../lib/SFML-2.5.1/include
*/

#include <SFML/Graphics.hpp>
#include "FPS.cpp"
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <omp.h>
#include <iostream>

//Representation of a particle, stores its physics information and size
struct Body
{
    sf::Vector2f pos;
    sf::Vector2f vel;
    sf::Vector2f acc;
    float mass = 1.0f; // 1kg
    float drag = 1.0f; // rho*C*Area � simplified drag for this example
	float radius = 5.0f; //size of each circle
    sf::CircleShape shape;

    Body() {
        this->pos = sf::Vector2f(0.0f, 0.0f);
        this->vel = sf::Vector2f(100.0f, 0.0f);
        this->acc = sf::Vector2f(0.0f, 0.0f);
        this->shape = sf::CircleShape(radius);
        this->shape.setOrigin(this->radius, this->radius);
        this->shape.setPosition(this->pos);
    }

    /**
     * Update pos and vel using "Velocity Verlet" integration
     * @param dt DeltaTime / time step [eg: 0.01]
     */
    void update(float dt)
    {
        sf::Vector2f new_pos = this->pos + this->vel * dt + this->acc * (dt * dt * 0.5f);
        sf::Vector2f new_acc = apply_forces(); // only needed if acceleration is not constant
        sf::Vector2f new_vel = this->vel + (this->acc + new_acc) * (dt * 0.5f);
		#pragma omp critical
        this->pos = new_pos;
		#pragma omp critical
        this->vel = new_vel;
        this->acc = new_acc;
        this->shape.setPosition(this->pos);
    }

    void setPosition(sf::Vector2f vec) 
    {
        this->pos = vec;
        this->shape.setPosition(this->pos);
    }

    sf::Vector2f apply_forces() const
    {
        sf::Vector2f grav_acc = sf::Vector2f( 0.0f, 100.0f); //Down is + in y, right is + in x
        sf::Vector2f drag_force = 0.5f * drag * vel; // D = 0.5 * (rho * C * Area * vel^2)
        sf::Vector2f drag_acc = drag_force / mass; // a = F/m
        return grav_acc- drag_acc;
    }
};

//Representation of a grid cell, contains a list of particles
struct Cell{
	std::vector<unsigned int> *cellParticles;
	int xPos;
	int yPos;
	int cellSize;

	Cell(int xPos, int yPos, int size){
		this->xPos = xPos;
		this->yPos = yPos;
		this->cellSize = size;
		this->cellParticles = new std::vector<unsigned int>();
	}

	void addParticle(unsigned int index){
		this->cellParticles->push_back(index);
	}

	void clearCell(){
		this->cellParticles->clear();
	}
};

static const int screenX = 800;
static const int screenY = 800;
static int cellSize; //should only ever be as small as the diameter of a particle
std::vector<std::vector<Cell*>> cells; //Matrix of all cells
std::vector<Body*> particles; //list of all particles

//function definitions
bool collide(Body* particle1, Body* particle2);
void solveCollision(Body* particle1, Body* particle2);
void findCollisions();
void check_cells_collisions(Cell* cell_1, Cell* cell_2);
void find_collisions_grid();
void updatePhysics(float dt);
void updatePhysicsSubtick(float dt, int subTicks);
void initializeCells();
Cell* getCell(sf::Vector2f);
void clearCells();
void fillCells();

int nThreads;

int main(int argc, char **argv)
{
	if (argc != 5)
    {
        std::cout << "Usage: export LD_LIBRARY_PATH=../lib/SFML-2.5.1/lib && ./app nThreads nParticles cellSize drawWindowFlag" << std::endl;
		std::cout << "cellSize must be divisible by 800 and >= 10, drawWindowFlag can be 0 or 1" << std::endl;
        exit(0);
    }
    
    nThreads = atoi(argv[1]);
    unsigned int nParticles = static_cast<unsigned int>(atoi(argv[2]));
	cellSize = atoi(argv[3]);
	bool drawWindow = atoi(argv[4]);

	cells = std::vector<std::vector<Cell*>>(screenX/cellSize, std::vector<Cell*>(screenY/cellSize));

    omp_set_num_threads(nThreads);

    sf::RenderWindow window(sf::VideoMode(800, 800), "SFML works!");
    window.setFramerateLimit(60);
    const float dt = 1.0f / static_cast<float>(60);
	window.setVisible(drawWindow);

	initializeCells();

    FPS fps;
    unsigned long long int counter = 0;
    int colorCounter = 0;
	int rgbCounter = 0;
	bool colorUp = true;

    sf::Font font; //set font
    if(!font.loadFromFile("./UbuntuMono-BI.ttf"))
    {
        std::cout << "Error loading font!" << std::endl;
    }

    sf::Text text;
    text.setFont(font);
    text.setCharacterSize(20);
    text.setFillColor(sf::Color::White);
    text.setPosition(150.f, 0.f);

    sf::Text text2;
    text2.setFont(font);
    text2.setCharacterSize(16);
    text2.setFillColor(sf::Color::White);
    text2.setPosition(150.f, 23.f);

    sf::Clock clock; //start clock
	sf::Time elapsed;

    while (window.isOpen()) //Main loop
    {        
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed){ //close the window if close is pressed
				window.close();
			}
                
        }
        if (counter%5 == 0 && particles.size() < nParticles) //controls rate of particle addition and changing colors
        {
            Body* newBody = new Body();
            if (colorCounter == 0) newBody->shape.setFillColor(sf::Color(255, 0, rgbCounter, 255));
            if (colorCounter == 1) newBody->shape.setFillColor(sf::Color(rgbCounter, 0, 255, 255));
            if (colorCounter == 2) newBody->shape.setFillColor(sf::Color(0, rgbCounter, 255, 255));
			if (colorCounter == 3) newBody->shape.setFillColor(sf::Color(0, 255, rgbCounter, 255));
            if (colorCounter == 4) newBody->shape.setFillColor(sf::Color(rgbCounter, 255, 0, 255));
            if (colorCounter == 5) newBody->shape.setFillColor(sf::Color(255, rgbCounter, 0, 255));
            if(colorUp)
			{
				if(rgbCounter<=254) rgbCounter+=6;
				else 
				{
					rgbCounter = 255;
					if(colorCounter == 5) colorCounter = 0;
					else colorCounter++;
					colorUp = false;
				}
			}
			else
			{
				if(rgbCounter>=0) rgbCounter-=25;
				else
				{
					rgbCounter = 0;
					if(colorCounter == 5) colorCounter = 0;
					else colorCounter++;
					colorUp = true;
				} 
			}
            particles.push_back(newBody);
        }

		clock.restart(); //restart clock
        updatePhysicsSubtick(dt, 8);
		elapsed += clock.getElapsedTime(); //measure time for update
		counter++;

        fps.update(); //FPS for display
        std::ostringstream ss;
        ss << fps.getFPS();

        window.setTitle(ss.str());

        window.clear(); //Draw all circles
        for (Body* particle : particles)
        {
            window.draw(particle->shape);
        }

        //display running info
		double avgElapsed = (elapsed.asMilliseconds()/(double)counter);
        text.setString("Avg ms/update: " + std::to_string(avgElapsed) + "    Number of Particles: " + std::to_string(particles.size()));
        window.draw(text);
        text2.setString("nThreads: " + std::to_string(nThreads));
        window.draw(text2);

        window.display();

        if(particles.size() == nParticles) //stop program when maximum number of particles is reached, comment out to keep running after
        {
            // std::cout << "FPS: " << fps.getFPS() << std::endl;
            window.close();
        }
    }

    //print results
    std::cout << "Total Elapsed time: " << elapsed.asSeconds() << " sec";
    std::cout << " Avg Time/physics update: " << (elapsed.asMilliseconds() / (double)counter) << " ms" << std::endl;

    particles.clear();
    return 0;
}

//Checks if two particles are colliding
//are colliding if distance between their centers is less than the sum of their radii
bool collide(Body* particle1, Body* particle2)
{
	const sf::Vector2f o2_o1 = particle1->pos - particle2->pos;
    const float dist2 = o2_o1.x * o2_o1.x + o2_o1.y * o2_o1.y;
    const float dist = sqrt(dist2);
    if(dist < (particle1->radius + particle2->radius))
        return true;
    else
        return false;
}

//Adjust two particles so that the no longer collide, does not fully separate particles
//as this makes simulation less smooth
void solveCollision(Body* particle1, Body* particle2)
{
    constexpr float response_coef = 1.0f; //Adjust how agressively particles are pushed apart, The higher the number the more jiggle but less squishing
    const sf::Vector2f o2_o1 = particle1->pos - particle2->pos;
    const float dist2 = o2_o1.x * o2_o1.x + o2_o1.y * o2_o1.y;
    
    const float dist = sqrt(dist2);
    const float delta = response_coef * 0.5f * ((particle1->radius+particle2->radius) - dist);
    //const float delta = dist / 4.0f;
    const sf::Vector2f col_vec = (o2_o1 / dist) * delta;
    sf::Vector2f newP1Pos = particle1->pos + col_vec;
    particle1->setPosition(newP1Pos);
    sf::Vector2f newP2Pos = particle2->pos - col_vec;
	particle2->setPosition(newP2Pos);
}

//Naive function for finding and solving all particle collisions
void findCollisions() 
{
	unsigned int i;
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

//Check collisions between all particles in two cells, 
//used to be used but now is called with same cell to do internal cell collisions
void check_cells_collisions(Cell* cell_1, Cell* cell_2)
{
     unsigned int i, j;
     auto cell1Particles = *cell_1->cellParticles;
     auto cell1ParticlesSize = cell1Particles.size();
     auto cell2Particles = *cell_2->cellParticles;
     auto cell2ParticlesSize = cell2Particles.size();
     //Essentially the naive approach but only per cell
     //split by threads
     #pragma omp parallel for 
     for (i = 0; i < cell1ParticlesSize;  i++)
     {
         for (j = 0; j < cell2ParticlesSize; j++)
         {
             auto particle1 = particles.at(cell1Particles[i]);
             auto particle2 = particles.at(cell2Particles[j]);
            
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

//Check collisions of every cell in the grid, split by threads
void find_collisions_grid()
{   
	unsigned int i, j;
    //loop through rows of grid
    #pragma omp parallel for private(j)
    for (i = 0; i < cells.size(); i++)
    {
        //loop through columns of grid
        for (j = 0; j < cells.at(i).size(); j++)
        {
            Cell* current_cell = cells.at(i).at(j);
            check_cells_collisions(current_cell, current_cell);
        }
    }
}

void updatePhysics(float dt)
{
    const float margin = 2.0f;
    unsigned int i;

	if(nThreads > 1){ //run grid approach if threads>1
		fillCells();
		find_collisions_grid();
	}

    #pragma omp parallel for
    for (i =0; i < particles.size(); i++)
    {
		Body* particle = particles.at(i);
		
        if(nThreads == 1){ //run naive approach if threads==1
			findCollisions();
		}
        
        if (particle->pos.x > screenX - margin - particle->radius) {
            particle->pos.x = screenX - margin - particle->radius;
            particle->setPosition(particle->pos);
			// particle->vel = sf::Vector2f(-particle->vel.x,particle->vel.y); //make particles bounce perfectly elasticly off window boundry
			
        }
        else if (particle->pos.x < margin + particle->radius) {
            particle->pos.x = margin + particle->radius;
            particle->setPosition(particle->pos);
			// particle->vel = sf::Vector2f(-particle->vel.x,particle->vel.y);
			
        }
        if (particle->pos.y > screenY - margin - particle->radius) {
            particle->pos.y = screenY - margin - particle->radius;
            particle->setPosition(particle->pos);
			// particle->vel = sf::Vector2f(particle->vel.x,-particle->vel.y);
			
        }
        else if (particle->pos.y < margin + particle->radius) {
            particle->pos.y = margin + particle->radius;
            particle->setPosition(particle->pos);
			// particle->vel = sf::Vector2f(particle->vel.x,-particle->vel.y);
			
        }
		particle->update(dt);
    }
}

//Splits one update into subTicks # of smaller updates
void updatePhysicsSubtick(float dt, int subTicks)
{
    const float sub_dt = dt / (float)subTicks;
    for (int i{ subTicks }; i--;)
    {
        updatePhysics(sub_dt);
    }
}

//Sets up grid of cells
void initializeCells(){
	unsigned int i, j;
	#pragma omp parallel for private(j)
	for(i = 0; i < cells.size(); i++)
	{
		int curX = screenX/cellSize * i;
		for(j = 0; j < cells.at(i).size(); j++)
		{
			int curY = screenY/cellSize * j;
			cells.at(i).at(j) = new Cell(curX, curY, cellSize);
		}
	}
}

//Gets the cell a position is inside
Cell* getCell(sf::Vector2f vec){
	unsigned int idx = static_cast<unsigned int>(vec.x/cellSize);
	unsigned int idy = static_cast<unsigned int>(vec.y/cellSize);
	return cells.at(idx).at(idy);
}

//Remove particles from all cells
void clearCells() {
    unsigned int i, j;
    #pragma omp parallel for private(j)
    for (i = 0; i < cells.size(); i++)
    {
        for (j = 0; j < cells.at(i).size(); j++)
        {
            cells.at(i).at(j)->clearCell();
        }
    }
}

//Fill cells with particles within their boundaries
void fillCells(){
	unsigned int i;
    clearCells();
    //seems consistently slower with parallel for
    // #pragma omp parallel for 
    for (i =0; i < particles.size(); i++)
	{
		Body* particle = particles.at(i);
		Cell* curCell = getCell(particle->pos);
		#pragma omp critical
		curCell->addParticle(i); //add all corners of particles to all cells it overlaps
        sf::Vector2f leftTop = sf::Vector2f(particle->pos.x - particle->radius, particle->pos.y - particle->radius);
        curCell = getCell(leftTop);
		#pragma omp critical
        curCell->addParticle(i);
        sf::Vector2f rightTop = sf::Vector2f(leftTop.x + particle->radius + particle->radius, leftTop.y);
        curCell = getCell(rightTop);
		#pragma omp critical
        curCell->addParticle(i);
        sf::Vector2f leftBottom = sf::Vector2f(particle->pos.x - particle->radius, particle->pos.y + particle->radius);
        curCell = getCell(leftBottom);
		#pragma omp critical
        curCell->addParticle(i);
        sf::Vector2f rightBottom = sf::Vector2f(leftBottom.x + particle->radius + particle->radius, leftBottom.y);
        curCell = getCell(rightBottom);
		#pragma omp critical
        curCell->addParticle(i);
	}
}