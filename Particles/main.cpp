/* TO COMPILE:
	g++ -fopenmp -g -c main.cpp -I../lib/SFML-2.5.1/include
	g++ main.o -fopenmp -o app -L../lib/SFML-2.5.1/lib -lsfml-graphics -lsfml-window -lsfml-system
   TO RUN:
	export LD_LIBRARY_PATH=../lib/SFML-2.5.1/lib && ./app

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

struct Body
{
    sf::Vector2f pos;
    sf::Vector2f vel; // 2 m/s along x-axis
    sf::Vector2f acc; // no acceleration at first
    float mass = 1.0f; // 1kg
    float drag = 1.0f; // rho*C*Area � simplified drag for this example
	float radius = 10.f;
    sf::CircleShape shape;

    Body() {
        this->pos = sf::Vector2f(0.0f, 0.0f);
        this->vel = sf::Vector2f(50.0f, 0.0f);
        this->acc = sf::Vector2f(0.0f, 0.0f);
        this->shape = sf::CircleShape(10.f);
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
        this->pos = new_pos;
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
        sf::Vector2f grav_acc = sf::Vector2f( 0.0f, 40.0f); //
        sf::Vector2f drag_force = 0.5f * drag * vel; // D = 0.5 * (rho * C * Area * vel^2)
        sf::Vector2f drag_acc = drag_force / mass; // a = F/m
        return grav_acc- drag_acc;
    }
};

struct Cell{
	static const unsigned int cellCapacity = 4;
	static const unsigned int maxCellIndex = cellCapacity-1;
	std::vector<unsigned int> *cellParticles;
	int xPos;
	int yPos;
	int cellSize;
	/*sf::FloatRect* rect;*/
    /*size_t particleCount = 0;
    size_t* cellParticles;*/

	Cell(int xPos, int yPos, int size){
		this->xPos = xPos;
		this->yPos = yPos;
		this->cellSize = size;
        //this->cellParticles = new size_t[cellCapacity];
		//this->rect = new sf::FloatRect(sf::Vector2f(xPos,yPos),sf::Vector2f(size,size));
		this->cellParticles = new std::vector<unsigned int>();
	}

	/*~Cell(){
		delete this->rect;
	}*/

	void addParticle(unsigned int index){
		 /*this->particleCount += particleCount < maxCellIndex;
         this->cellParticles[particleCount] = index;*/
		this->cellParticles->push_back(index);
	}

	void clearCell(){
		this->cellParticles->clear();
        /*this->particleCount = 0;*/
	}
};

static const int screenX = 800;
static const int screenY = 800;
static const int cellSize = 40;
std::vector<std::vector<Cell*>> cells(screenX/cellSize, std::vector<Cell*>(screenY/cellSize));
std::vector<Body*> particles;

bool collide(Body* particle1, Body* particle2);
void solveCollision(Body* particle1, Body* particle2);
void findCollisions();
void findCollisionsInCell(Cell* cell);
void check_cells_collisions(Cell* cell_1, Cell* cell_2);
void find_collisions_grid();
void updatePhysics(float dt);
void updatePhysicsSubtick(float dt, int subTicks);
void initializeCells();
Cell* getCell(float xPos, float yPos);
void clearCells();
void fillCells();

int main(int argc, char **argv)
{
	if (argc != 3)
    {
        std::cout << "Usage: export LD_LIBRARY_PATH=../lib/SFML-2.5.1/lib && ./app nThreads nParticles" << std::endl;
        exit(0);
    }
    
    int nThreads = atoi(argv[1]);
    int nParticles = atoi(argv[2]);

    omp_set_num_threads(nThreads);

    sf::RenderWindow window(sf::VideoMode(800, 800), "SFML works!");
    window.setFramerateLimit(60);
    const float dt = 1.0f / static_cast<float>(60);

    particles.push_back(new Body());
	//std::cout << particles.at(0)->radius << std::endl;

	initializeCells();
    fillCells();

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

    sf::Text text2;
    text2.setFont(font);
    text2.setCharacterSize(16);
    text2.setFillColor(sf::Color::White);
    text2.setPosition(150.f, 23.f);

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
        if (counter % 15 == 0 && particles.size() < nParticles)
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
        //updatePhysics(dt);

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
        text2.setString("nThreads: " + std::to_string(nThreads));
        window.draw(text2);

        window.display();
        if (counter < 60) counter++;
        else counter = 0;
    }

    sf::Time elapsed = clock.getElapsedTime();
    std::cout << "Elapsed time: " << elapsed.asSeconds() << std::endl;

    particles.clear();
    return 0;
}

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

void solveCollision(Body* particle1, Body* particle2)
{
    constexpr float response_coef = 1.0f;
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

void check_cells_collisions(Cell* cell_1, Cell* cell_2)
{
     auto cell1Particles = *cell_1->cellParticles;
     auto cell1ParticlesSize = cell1Particles.size();
     auto cell2Particles = *cell_2->cellParticles;
     auto cell2ParticlesSize = cell2Particles.size();
    //  #pragma omp parallel for 
     for (unsigned int i = 0; i < cell1ParticlesSize;  i++)
     {
         for (unsigned int j = 0; j < cell2ParticlesSize; j++)
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

void find_collisions_grid()
{   
	unsigned int i, j;
    //loop through rows of grid (skip top and bottom rows)
    // #pragma omp parallel for private(j)
    for (i = 0; i < cells.size(); i++)
    {
        //loop through columns of grid (skip left and right columns)
        for (j = 0; j < cells.at(i).size(); j++)
        {
            Cell* current_cell = cells.at(i).at(j);
            //Iterate on all surrounding cells, including current one
            if (i == 0 || j == 0 || i == cells.size() - 1 || j == cells.at(i).size() - 1) 
            {
                check_cells_collisions(current_cell, current_cell);
            }
            else
            {
                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        Cell* other_cell = cells.at(i + dx).at(j + dy);
                        check_cells_collisions(current_cell, other_cell);
                    }
                }
            }
        }
    }
}

void updatePhysics(float dt)
{
    const float margin = 2.0f;
    unsigned int i;

	fillCells();
	find_collisions_grid();

    #pragma omp parallel for
    for (i =0; i < particles.size(); i++)
    {
		Body* particle = particles.at(i);
		
        // findCollisions();
        
        if (particle->pos.x > screenX - margin - particle->radius) {
            particle->pos.x = screenX - margin - particle->radius;
            particle->setPosition(particle->pos); //remove for expanseVer
			particle->vel = sf::Vector2f(-particle->vel.x,particle->vel.y);
			
        }
        else if (particle->pos.x < margin + particle->radius) {
            particle->pos.x = margin + particle->radius;
            particle->setPosition(particle->pos); //remove for expanseVer
			particle->vel = sf::Vector2f(-particle->vel.x,particle->vel.y);
			
        }
        if (particle->pos.y > screenY - margin - particle->radius) {
            particle->pos.y = screenY - margin - particle->radius;
            particle->setPosition(particle->pos); //remove for expanseVer
			particle->vel = sf::Vector2f(particle->vel.x,-particle->vel.y);
			
        }
        else if (particle->pos.y < margin + particle->radius) {
            particle->pos.y = margin + particle->radius;
            particle->setPosition(particle->pos); //remove for expanseVer
			particle->vel = sf::Vector2f(particle->vel.x,-particle->vel.y);
			
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

Cell* getCell(float xPos, float yPos){
	unsigned int idx = static_cast<unsigned int>(xPos/cellSize);
	unsigned int idy = static_cast<unsigned int>(yPos/cellSize);
	return cells.at(idx).at(idy);
}

// std::vector<Cell*> getSurroundingCells(Cell* centerCell){
// 	auto posx = centerCell->xPos;
// 	auto posy = centerCell->yPos;
// 	int idx = std::floor(posx/cellSize);
// 	int idy = std::floor(posy/cellSize);
// 	std::vector<Cell*> result;
// 	if(idx - 1 >=0 && idy -1 >=0) result.push_back(cells.at(idx - 1).at(idy -1)); //top left
// 	if(idy -1 >=0) result.push_back(cells.at(idx).at(idy -1)); //top center
// 	if(idx +1<=cells.size() && idy -1 >=0) result.push_back(cells.at(idx+1).at(idy -1)); //top right
// 	if(idx-1 >=0) result.push_back(cells.at(idx-1).at(idy)); //left
// 	if(idx+1 <=cells.size()) result.push_back(cells.at(idx+1).at(idy)); //right
// 	if(idx - 1 >=0 && idy +1 <=cells.at(idx - 1).size()) result.push_back(cells.at(idx - 1).at(idy +1)); //top left
// 	if(idy +1 <=cells.at(idx).size()) result.push_back(cells.at(idx).at(idy +1)); //top center
// 	if(idx +1<=cells.size() && idy +1 <=cells.at(idx + 1).size()) result.push_back(cells.at(idx+1).at(idy +1)); //top right
// 	return result;
// }

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

void fillCells(){
	unsigned int i;
    clearCells();
    //#pragma omp parallel for
    for (i =0; i < particles.size(); i++)
	{
		Body* particle = particles.at(i);
		Cell* curCell = getCell(particle->pos.x, particle->pos.y);
		curCell->addParticle(i);
        //auto particlePos = particle->shape.getGlobalBounds();
		/*auto surroundingCells = getSurroundingCells(curCell);
		for(auto cell : surroundingCells)
		{
			auto cellRect = *cell->rect;
			if(particlePos.intersects(cellRect)){
				cell->addParticle(i);
			}
		}*/
		
	}
}