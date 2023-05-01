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

struct Cell{
	// static const size_t cellCapacity = 4;
	// static const size_t maxCellIndex = cellCapacity-1;
	std::vector<size_t> *cellParticles;
	size_t xPos;
	size_t yPos;
	int cellSize;
	sf::FloatRect* rect;

	Cell(size_t xPos, size_t yPos, int size){
		this->xPos = xPos;
		this->yPos = yPos;
		this->cellSize = size;
		this->rect = new sf::FloatRect(sf::Vector2f(xPos,yPos),sf::Vector2f(size,size));
		this->cellParticles = new std::vector<size_t>();
	}

	~Cell(){
		delete this->rect;
	}

	// size_t particleCount = 0;
	// size_t particles[cellCapacity];

	void addParticle(size_t index){
		// particleCount += particleCount < maxCellIndex;
		// particles[particleCount] = index;
		this->cellParticles->push_back(index);
	}

	void clearCell(){
		this->cellParticles->clear();
	}
}; //cells are 80px by 80px --resize by window and particle size

static const size_t screenX = 800;
static const size_t screenY = 800;
static const size_t cellSize = 80;
std::vector<std::vector<Cell*>> cells(screenX/cellSize, std::vector<Cell*>(screenY/cellSize));
std::vector<Body*> particles;
void findCollisions();
void updatePhysicsSubtick(float dt, int subTicks);
void initializeCells();
void fillCells();
Cell* getCell(float xPos, float yPos);
std::vector<Cell*> getSurroundingCells(Cell* centerCell);

int main()
{
    sf::RenderWindow window(sf::VideoMode(800, 800), "SFML works!");
    window.setFramerateLimit(60);
    const float dt = 1.0f / static_cast<float>(60);
    /*sf::CircleShape shape(10.f);
    shape.setFillColor(sf::Color::Green);*/
    particles.push_back(new Body());
	std::cout << particles.at(0)->shape.getRadius();
	initializeCells();
	fillCells();
    /*particles.push_back(new Body());
    particles.at(1)->setPosition(400.0f, 0.0f);
    particles.at(1)->vel.x = -40.0f;
    particles.at(1)->shape.setFillColor(sf::Color(0,255,0,255));*/
    //particle1.setPosition(400.0f,400.0f);

    FPS fps;
    unsigned int counter = 0;
    int colorCounter = 0;
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed){
				window.close();
			}
                
        }
        if (counter % 15 == 0 && particles.size() < 35)
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
        window.display();
        if (counter < 60) counter++;
        else counter = 0;
    }

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

//void solveCollision(Body* particle1, Body* particle2)
//{
//    float overlap = (particle1->shape.getRadius() + particle2->shape.getRadius()) - (float)distance(particle1->pos.x, particle1->pos.y, particle2->pos.x, particle2->pos.y);
//    sf::Vector2f vec = particle1->shape.getPosition() - particle2->shape.getPosition(); //particle 2 is start, particle 1 is terminal direction is toward particle 1
//    vec /= overlap;
//    vec *= overlap/2;
//    particle2->shape.move(-vec);
//    particle2->updatePositionByShape();
//    particle1->shape.move(vec);
//    particle1->updatePositionByShape();
//}

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

void findCollisionsCell(Cell* cell) 
{
	int i, j;
	#pragma omp parallel for private(j)
    for (i =0; i < cell->cellParticles->size(); i++)
    {
		auto particle1 = particles.at(cell->cellParticles->at(i));
		for (j =0; j < cell->cellParticles->size(); i++)
		{
			auto particle2 = particles.at(cell->cellParticles->at(j));
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

void findCollisionsCells() 
{
	int i, j;
	#pragma omp parallel for private(j)
	for(i = 0; i < cells.size(); i++)
	{
		for(j = 0; j < cells.at(i).size(); j++)
		{
			auto cell = cells.at(i).at(j);
			findCollisionsCell(cell);
		}
	}
}

void check_cells_collisions(Cell* cell_1, Cell* cell_2)
{
    // #pragma omp parallel for
    // for (auto& particle1Idx : *cell_1->particles)
    // {
    //     for (auto& particle2Idx : *cell_2->particles)
    //     {
    //         auto particle1 = particles.at(particle1Idx);
    //         auto particle2 = particles.at(particle2Idx);
            
    //         if (particle1 != particle2)
    //         {
    //             if (collide(particle1, particle2))
    //             {
    //                 solveCollision(particle1, particle2);
    //             }
    //         }
    //     }
    // }
}

void find_collisions_grid()
{   
	int i, j;
    //loop through rows of grid (skip top and bottom rows)
    // #pragma omp parallel for private(j)
    for (i = 1; i < cells.size() - 1; i++)
    {
        //loop through columns of grid (skip left and right columns)
        for (j = 1; j < cells.at(i).size() - 1; j++)
        {
            Cell* current_cell = getCell(i, j);
            //Iterate on all surrounding cells, including current one
            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    Cell* other_cell = getCell(i + dx, j + dy);
                    check_cells_collisions(current_cell, other_cell);
                }
            }
        }
    }
}

void updatePhysics(float dt)
{
    const float margin = 2.0f;
    int i;

	fillCells();
	// find_collisions_grid();
	// findCollisionsCells();

    #pragma omp parallel for
    for (i =0; i < particles.size(); i++)
    {
		auto particle = particles.at(i);
		
        // findCollisions();
        
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

void initializeCells(){
	int i, j;
	#pragma omp parallel for private(j)
	for(i = 0; i < cells.size(); i++)
	{
		size_t curX = screenX/cellSize * i;
		for(j = 0; j < cells.at(i).size(); j++)
		{
			size_t curY = screenY/cellSize * j;
			cells.at(i).at(j) = new Cell(curX, curY, cellSize);
		}
	}
}

Cell* getCell(float xPos, float yPos){
	int idx = std::floor(xPos/cellSize);
	int idy = std::floor(yPos/cellSize);
	return cells.at(idx).at(idy);
}

std::vector<Cell*> getSurroundingCells(Cell* centerCell){ //BUGGED
	auto posx = centerCell->xPos;
	auto posy = centerCell->yPos;
	int idx = std::floor(posx/cellSize);
	int idy = std::floor(posy/cellSize);
	std::vector<Cell*> result;
	if(idx - 1 >=0 && idy -1 >=0) result.push_back(cells.at(idx - 1).at(idy -1)); //top left
	if(idy -1 >=0) result.push_back(cells.at(idx).at(idy -1)); //top center
	if(idx +1<=cells.size() && idy -1 >=0) result.push_back(cells.at(idx+1).at(idy -1)); //top right
	if(idx-1 >=0) result.push_back(cells.at(idx-1).at(idy)); //left
	if(idx+1 <=cells.size()) result.push_back(cells.at(idx+1).at(idy)); //right
	if(idx - 1 >=0 && idy +1 <=cells.at(idx - 1).size()) result.push_back(cells.at(idx - 1).at(idy +1)); //top left
	if(idy +1 <=cells.at(idx).size()) result.push_back(cells.at(idx).at(idy +1)); //top center
	if(idx +1<=cells.size() && idy +1 <=cells.at(idx + 1).size()) result.push_back(cells.at(idx+1).at(idy +1)); //top right
	return result;
}

void fillCells(){
	int i, j;

    // #pragma omp parallel for
    for (i =0; i < particles.size(); i++)
	{
		auto particle = particles.at(i);
		auto particlePos = particle->shape.getGlobalBounds();
		auto curCell = getCell(particle->pos.x, particle->pos.y);
		curCell->addParticle(i);
		auto surroundingCells = getSurroundingCells(curCell);
		for(auto cell : surroundingCells)
		{
			auto cellRect = *cell->rect;
			if(particlePos.intersects(cellRect)){
				cell->addParticle(i);
			}
		}
		
	}
}