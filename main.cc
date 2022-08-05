#include <SDL/SDL.h>
#include <math.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <random>
#include <stack>
#include <queue>
#include <cassert>
#include <sys/time.h>

struct rand_context
{
private:
	std::mt19937 mt;
public:
	rand_context()
	{
		mt = std::mt19937(std::random_device()());
	}
	float randf(float min, float max)
	{
		return std::uniform_real_distribution<float>(min, max)(mt);
	}
};

union v2f
{
	struct { float a, b; };
	struct { float x, y; };
	struct { float v[2]; };

	inline v2f operator+(const v2f &b) const
	{
		return {x + b.x, y + b.y};
	}
	inline v2f &operator+=(const v2f &b)
	{
		x += b.x;
		y += b.y;
		return *this;
	}
	inline v2f operator-(const v2f &b) const
	{
		return {x - b.x, y - b.y};
	}
	inline v2f &operator-=(const v2f &b)
	{
		x -= b.x;
		y -= b.y;
		return *this;
	}
	inline v2f operator*(const v2f &b) const
	{
		return {x * b.x, y * b.y};
	}
	inline v2f operator*(const float b) const
	{
		return {x * b, y * b};
	}
	inline v2f &operator*=(const v2f &b)
	{
		x *= b.x;
		y *= b.y;
		return *this;
	}
	inline v2f &operator*=(const float b)
	{
		x *= b;
		y *= b;
		return *this;
	}
	inline v2f operator/(const v2f &b) const
	{
		return {x / b.x, y / b.y};
	}
	inline v2f operator/(const float b) const
	{
		return {x / b, y / b};
	}
	inline v2f &operator/=(const v2f &b)
	{
		x /= b.x;
		y /= b.y;
		return *this;
	}
	inline v2f &operator/=(const float b)
	{
		x /= b;
		y /= b;
		return *this;
	}
	inline v2f operator~() const
	{
		float il = 1.0f / std::sqrt(x*x + y*y);
		return {x * il, y * il};
	}
	inline float dot(const v2f &b) const
	{
		return x*b.x + y*b.y;
	}
	inline float length() const
	{
		return std::sqrt(x*x + y*y);
	}
	inline float length2() const
	{
		return x*x + y*y;
	}
	inline v2f mask_xy() const
	{
		return *this;
	}
	inline v2f mask_yx() const
	{
		return {y, x};
	}
	inline v2f mask_x0() const
	{
		return {x, 0.0f};
	}
	inline v2f mask_y0() const
	{
		return {y, 0.0f};
	}
	inline v2f mask_0x() const
	{
		return {0.0f, x};
	}
	inline v2f mask_0y() const
	{
		return {0.0f, y};
	}
};

struct m3f
{
	float m[3][3]{ 0 };
	
	static inline m3f identity()
	{
		return {{
			1.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 1.0f
		}};
	}
	static inline m3f scale(float x, float y)
	{
		return {{
			   x, 0.0f, 0.0f,
			0.0f,    y, 0.0f,
			0.0f, 0.0f, 1.0f
		}};
	}
	static inline m3f translate(float x, float y)
	{
		return {{
			1.0f, 0.0f,    x,
			0.0f, 1.0f,    y,
			0.0f, 0.0f, 1.0f
		}};
	}
	inline m3f operator*(const m3f &b)
	{
		m3f o;
		for (uint64_t i = 0; i < 3; i++)
			for (uint64_t j = 0; j < 3; j++)
				o.m[i][j] = m[i][0] * b.m[0][j] + m[i][1] * b.m[1][j] + m[i][2] * b.m[2][j];
		return o;
	}
	inline v2f operator*(const v2f &b)
	{
		v2f o;
		o.x = m[0][0] * b.x + m[0][1] * b.y + m[0][2];
		o.y = m[1][0] * b.x + m[1][1] * b.y + m[1][2];
		return o;
	}
	inline m3f fast_inverse()
	{
		/*
		 * [a, 0, c]^-1   [1/a,   0, -c/a]
		 * [0, b, d]    = [  0, 1/b, -d/b]
		 * [0, 0, 1]      [  0,   0,    1]
		 */
		float ia = 1.0f / m[0][0];
		float ib = 1.0f / m[1][1];
		return {{
			  ia, 0.0f, -m[0][2] * ia,
			0.0f,   ib, -m[1][2] * ib,
			0.0f, 0.0f,          1.0f
		}};
	}
};

struct rgb
{
	float r, g, b;
	inline void to_3u8(uint8_t &ro, uint8_t &go, uint8_t &bo)
	{
		ro = (uint8_t)std::max(0.0f, std::min(r * 256.0f, 255.0f));
		go = (uint8_t)std::max(0.0f, std::min(g * 256.0f, 255.0f));
		bo = (uint8_t)std::max(0.0f, std::min(b * 256.0f, 255.0f));
	}
};

template <uint16_t W, uint16_t H>
struct screen
{
	rgb pixels[W*H];
	screen()
	{
		for (uint64_t i = 0; i < W*H; i++)
		{
			pixels[i].r = 0.0f;
			pixels[i].g = 0.0f;
			pixels[i].b = 0.0f;
		}
	}
	inline uint16_t width()
	{
		return W;
	}
	inline uint16_t height()
	{
		return H;
	}
	inline rgb get(int x, int y)
	{
		if (x < 0 || x >= W || y < 0 || y >= H) return {0.0f, 0.0f, 0.0f};
		return pixels[y * W + x];
	}
	inline void set(int x, int y, rgb col)
	{
		if (x < 0 || x >= W || y < 0 || y >= H) return;
		pixels[y * W + x] = col;
	}
	inline void set_scanline(int x1, int x2, int y, rgb col)
	{
		if (y < 0 || y >= H) return;
		x1 = std::max(x1, 0);
		x2 = std::min(x2, W - 1);
		for (int i = x1; i <= x2; i++)
			pixels[y * W + i] = col;
	}
	inline void fill(rgb col)
	{
		for (uint64_t i = 0; i < W * H; i++)
			pixels[i] = col;
	}
};

template <uint16_t W, uint16_t H>
void draw_circle(screen<W, H> *scr, int cx, int cy, float r, rgb col)
{
	int x = 0;
	int y = r;
	while (x <= y)
	{
		scr->set(cx + x, cy + y, col);
		scr->set(cx + y, cy + x, col);
		scr->set(cx + x, cy - y, col);
		scr->set(cx + y, cy - x, col);
		scr->set(cx - x, cy + y, col);
		scr->set(cx - y, cy + x, col);
		scr->set(cx - x, cy - y, col);
		scr->set(cx - y, cy - x, col);
		x++;
		if ((float)(x*x + y*y) > r*r) y--;
	}
}

template <uint16_t W, uint16_t H>
void fill_circle(screen<W, H> *scr, int cx, int cy, float r, rgb col)
{
	int x = r;
	int y = 0;
	while (x >= y)
	{
		scr->set_scanline(cx - x, cx + x, cy - y, col);
		scr->set_scanline(cx - x, cx + x, cy + y, col);
		y++;
		if ((float)(x*x + y*y) > r*r)
		{
			scr->set_scanline(cx - (y - 1), cx + (y - 1), cy + x, col);
			scr->set_scanline(cx - (y - 1), cx + (y - 1), cy - x, col);
			x--;
		}
	}
}

template <uint16_t W, uint16_t H>
void draw_rect(screen<W, H> *scr, int x1, int y1, int x2, int y2, rgb col)
{
	if (y1+1 > y2) std::swap(y1, y2);
	if (x1 > x2) std::swap(x1, x2);
	scr->set_scanline(x1, x2, y1, col);
	scr->set_scanline(x1, x2, y2, col);
	for (int i = y1+1; i < y2; i++)
	{
		scr->set(x1, i, col);
		scr->set(x2, i, col);
	}
}

template <uint16_t W, uint16_t H>
void draw_line(screen<W, H> *scr, int x1, int y1, int x2, int y2, rgb col)
{
	float dx = x2 - x1;
	float dy = y2 - y1;
	int s = std::max(std::fabs(dx), std::fabs(dy));
	dx /= s;
	dy /= s;
	for (int i = 0; i < s; i++)
		scr->set(x1 + dx*i, y1 + dy*i, col);
}

struct rect
{
	v2f pos, size;
	inline bool contains(v2f p) const
	{
		return p.x >= pos.x && p.x < pos.x + size.x && p.y >= pos.y && p.y < pos.y + size.y;
	}
	inline bool contains(rect r) const
	{
		return r.pos.x >= pos.x && r.pos.x + r.size.x < pos.x + size.x && r.pos.y >= pos.y && r.pos.y + r.size.y < pos.y + size.y;
	}
	inline bool intersects(rect r) const
	{
		return r.pos.x + r.size.x >= pos.x && r.pos.x < pos.x + size.x && r.pos.y + r.size.y >= pos.y && r.pos.y < pos.y + size.y;
	}
};

template <typename T>
struct quad_tree_node;

template <typename T>
struct obj_with_area
{
	using list_iter = std::_List_iterator<obj_with_area<T> *>;
	T o;
	rect r;
	quad_tree_node<T> *p;
	list_iter it;
};

template<typename T>
struct quad_tree_node
{
	std::list<obj_with_area<T> *> l;
	quad_tree_node<T> *n[4] = { 0 };
	quad_tree_node<T> *p = 0;
	rect node_rect;
	uint64_t len;
	uint64_t depth;
};

template <typename T, int DEPTH> requires requires(T o) { o.get_bounding_box(); }
class quad_tree
{
public:
	using list_iter = std::_List_iterator<obj_with_area<T> *>;

	quad_tree_node<T> *root;
	quad_tree(const rect &world_rect)
	{
		this->root = new quad_tree_node<T>;
		this->root->node_rect = world_rect;
		this->root->len = 0;
		this->root->depth = 0;
	}
	static void clear_child(quad_tree_node<T> *p)
	{
		for (uint64_t i = 0; i < 4; i++)
		{
			if (p->n[i] != nullptr)
			{
				quad_tree::clear_child(p->n[i]);
				for (auto q : p->n[i]->l) delete q;
				p->n[i]->l.clear();
				delete p->n[i];
				p->n[i] = nullptr;
			}
		}
	}
	~quad_tree()
	{
		quad_tree::clear_child(root);
		for (auto p : root->l) delete p;
		delete root;
	}
	obj_with_area<T> *insert(T &&o)
	{
		rect r = o.get_bounding_box();
		assert(root->node_rect.contains(r));

		quad_tree_node<T> *parent = root;
		rect child_rect[4];
		int depth = 0;
		while (depth < DEPTH)
		{
			v2f hs = parent->node_rect.size * 0.5f;
			child_rect[0].pos = parent->node_rect.pos;
			child_rect[0].size = hs;
			child_rect[1].pos = parent->node_rect.pos + hs.mask_x0();
			child_rect[1].size = hs;
			child_rect[2].pos = parent->node_rect.pos + hs.mask_0y();
			child_rect[2].size = hs;
			child_rect[3].pos = parent->node_rect.pos + hs;
			child_rect[3].size = hs;
			int next;
			int num_contains = 0;
			for (uint64_t i = 0; i < 4; i++)
			{
				if (child_rect[i].contains(r))
				{
					next = i;
					num_contains++;
				}
			}
			if (num_contains == 1)
			{
				if (parent->n[next])
					parent = parent->n[next];
				else
				{
					for (uint64_t i = 0; i < 4; i++)
					{
						parent->n[i] = new quad_tree_node<T>;
						parent->n[i]->node_rect = child_rect[i];
						parent->n[i]->p = parent;
						parent->n[i]->len = 0;
						parent->n[i]->depth = parent->depth + 1;
					}
					parent = parent->n[next];
				}
			}
			else break;
			depth++;
		}
		for (quad_tree_node<T> *q = parent; q != nullptr; q = q->p)
			q->len++;
		parent->l.push_back(new obj_with_area<T>{o, r, parent});
		parent->l.back()->it = --parent->l.end();
		return parent->l.back();
	}
	std::list<obj_with_area<T> *> get_intersect(const rect &r)
	{
		std::list<obj_with_area<T> *> out;
		std::stack<quad_tree_node<T> *> to_process;
		to_process.push(root);
		while (!to_process.empty())
		{
			quad_tree_node<T> *p = to_process.top();
			to_process.pop();
			if (r.intersects(p->node_rect))
			{
				for (uint64_t i = 0; i < 4; i++)
					if (p->n[i] != nullptr)
						to_process.push(p->n[i]);
				for (auto o : p->l)
				{
					if (r.intersects(o->r))
						out.push_back(o);
				}
			}
		}
		return out;
	}
	std::list<obj_with_area<T> *> get_contains(const v2f &v)
	{
		std::list<obj_with_area<T> *> out;
		std::stack<quad_tree_node<T> *> to_process;
		to_process.push(root);
		while (!to_process.empty())
		{
			quad_tree_node<T> *p = to_process.top();
			to_process.pop();
			if (p->node_rect.contains(v))
			{
				for (uint64_t i = 0; i < 4; i++)
					if (p->n[i] != nullptr)
						to_process.push(p->n[i]);
				for (auto o : p->l)
				{
					if (o->r.contains(v))
						out.push_back(o);
				}
			}
		}
		return out;
	}
	void detach_obj(obj_with_area<T> *o)
	{
		quad_tree_node<T> *p = o->p;
		p->l.erase(o->it);
		o->p = nullptr;
		quad_tree_node<T> *q = p;
		while (q != nullptr)
		{
			q->len--;
			q = q->p;
		}
		q = p;
		while (q != nullptr)
		{
			if (q->l.size() == q->len)
			{
				quad_tree::clear_child(q);
			}
			q = q->p;
		}
	}
	void reattach_obj(obj_with_area<T> *o)
	{
		assert(root->node_rect.contains(o->r));
		quad_tree_node<T> *parent = root;
		rect child_rect[4];
		int depth = 0;
		while (depth < DEPTH)
		{
			v2f hs = parent->node_rect.size * 0.5f;
			child_rect[0].pos = parent->node_rect.pos;
			child_rect[0].size = hs;
			child_rect[1].pos = parent->node_rect.pos + hs.mask_x0();
			child_rect[1].size = hs;
			child_rect[2].pos = parent->node_rect.pos + hs.mask_0y();
			child_rect[2].size = hs;
			child_rect[3].pos = parent->node_rect.pos + hs;
			child_rect[3].size = hs;
			int next;
			int num_contains = 0;
			for (uint64_t i = 0; i < 4; i++)
			{
				if (child_rect[i].contains(o->r))
				{
					next = i;
					num_contains++;
				}
			}
			if (num_contains == 1)
			{
				if (parent->n[next])
					parent = parent->n[next];
				else
				{
					for (uint64_t i = 0; i < 4; i++)
					{
						parent->n[i] = new quad_tree_node<T>;
						parent->n[i]->node_rect = child_rect[i];
						parent->n[i]->p = parent;
						parent->n[i]->len = 0;
						parent->n[i]->depth = parent->depth + 1;
					}
					parent = parent->n[next];
				}
			}
			else break;
			depth++;
		}
		for (quad_tree_node<T> *q = parent; q != nullptr; q = q->p)
			q->len++;
		parent->l.push_back(o);
		parent->l.back()->p = parent;
		parent->l.back()->it = --parent->l.end();
	}
	void update_obj(obj_with_area<T> *o)
	{
		detach_obj(o);
		o->r = o->o.get_bounding_box();
		return reattach_obj(o);
	}
};

struct ball
{
	v2f p, v, a;
	float r;
	bool is_moving;
	ball(v2f p, float r)
	{
		this->p = p;
		this->v = {0.0f, 0.0f};
		this->a = {0.0f, 0.0f};
		this->r = r;
		this->is_moving = false;
	}
	rect get_bounding_box() const
	{
		rect o;
		v2f hs = {r, r};
		o.pos = p - hs;
		o.size = hs * 2.0f;
		return o;
	}
};

void alt_main()
{
	v2f p = {0.0f, 0.0f}, s = {32.0f, 32.0f};
	rect r = {p, s};
	quad_tree<ball, 5> *qt = new quad_tree<ball, 5>(r);
	int c;
	for (;;)
	{
		std::cout << "1. create new ball\n";
		std::cout << "2. get intersect\n";
		std::cout << "3. update\n";
		std::cout << "4. trap\n";
		std::cout << "5. exit" << std::endl;
		std::cin >> c;
		switch (c)
		{
		case 1:
		{
			std::cout << "ball pos: ";
			v2f p;
			std::cin >> p.x >> p.y;
			std::cout << "ball radius: ";
			float r;
			std::cin >> r;
			qt->insert({p, r});
			break;
		}
		case 2:
		{
			std::cout << "rect pos: ";
			v2f p;
			std::cin >> p.x >> p.y;
			std::cout << "rect size: ";
			v2f s;
			std::cin >> s.x >> s.y;
			rect r = {p, s};
			auto l = qt->get_intersect(r);
			for (auto o : l)
			{
				std::cout << "pos: (" << o->o.p.x << ", " << o->o.p.y << ") radius: " << o->o.r << '\n';
			}
			break;
		}
		case 3:
		{
			std::cout << "rect pos: ";
			v2f p;
			std::cin >> p.x >> p.y;
			std::cout << "rect size: ";
			v2f s;
			std::cin >> s.x >> s.y;
			rect r = {p, s};
			auto l = qt->get_intersect(r);
			if (l.empty())
				std::cout << "no balls intersect the rectangle" << std::endl;
			else
			{
				auto &b = l.front()->o;
				std::cout << "new pos: ";
				v2f p;
				std::cin >> p.x >> p.y;
				std::cout << "new radius: ";
				float r;
				std::cin >> r;
				std::cout << "ball (" << b.p.x << ", " << b.p.y << ", " << b.r << ")";
				b.p = p;
				b.r = r;
				qt->update_obj(l.front());
				std::cout << " -> (" << b.p.x << ", " << b.p.y << ", " << b.r << ")" << std::endl;
			}
			break;
		}
		case 4: { asm("int3"); break; }
		case 5: { exit(0); break; }
		default:
			break;
		}
	}
	delete qt;
}


int main()
{
	constexpr uint16_t W = 1280, H = 720;
	constexpr float as = (float)W / H;
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Surface *surface = SDL_SetVideoMode(W, H, 32, 0);
	screen<W, H> *scr = new screen<W, H>;
	rand_context rc;
	uint8_t *keys = SDL_GetKeyState(nullptr);
	uint8_t prev_mouse_state = 0;

	const float world_size = 1024.0f;
	const rect world_rect = {{-world_size, -world_size}, {world_size * 2, world_size * 2}};
	constexpr int MAX_DEPTH = 12;
	quad_tree<ball, MAX_DEPTH> world(world_rect);
	for (uint64_t i = 0; i < 30000; i++)
	{
		ball b({rc.randf(-1000, 1000), rc.randf(-1000, 1000)}, rc.randf(0.25f, 4.0f));
		b.is_moving = false;
		world.insert(std::move(b));
	}
	obj_with_area<ball> *ball_move = nullptr;
	bool ball_moving = false;

	float zoom = 0.125f;
	v2f zoom_center = {0.0f, 0.0f};

	struct timeval t1, t2;
	gettimeofday(&t1, NULL);

	bool is_running = true;
	for (uint64_t t = 0; is_running; t++)
	{
		gettimeofday(&t2, NULL);
		int64_t sec = t2.tv_sec - t1.tv_sec;
		int64_t usec = t2.tv_usec - t1.tv_usec;
		t1 = t2;
		float dt = (float)sec + (float)usec / 1000000.0f;
		char s[100];
		sprintf(s, "fps: %.3f", 1.0f / dt);
		SDL_WM_SetCaption(s, NULL);

		for (SDL_Event e; SDL_PollEvent(&e);)
		{
			if (e.type == SDL_QUIT) is_running = 0;
		}
		int mouse_x, mouse_y;
		uint8_t mouse_state = SDL_GetMouseState(&mouse_x, &mouse_y);

		v2f move = {0.0f, 0.0f};
		const v2f right = {1.0f, 0.0f};
		const v2f up = {0.0f, 1.0f};
		if (keys[SDLK_LEFT]) move -= right;
		if (keys[SDLK_RIGHT]) move += right;
		if (keys[SDLK_UP]) move += up;
		if (keys[SDLK_DOWN]) move -= up;
		if (move.length2() != 0.0f)
		{
			zoom_center += ~move * dt / zoom;
		}
		if (keys[SDLK_z]) zoom *= std::pow(2.0f, dt);
		if (keys[SDLK_x]) zoom *= std::pow(0.5f, dt);

		scr->fill({0.0f, 0.0f, 0.0f});

		m3f camera_mat = m3f::translate(zoom_center.x, zoom_center.y) * m3f::scale(1.0f / zoom, 1.0f / zoom);
		v2f camera_bl = {-as*0.5f, -0.5f};
		v2f camera_tr = camera_bl * -1.0f;
		camera_bl = camera_mat * camera_bl;
		camera_tr = camera_mat * camera_tr;
		v2f camera_size = camera_tr - camera_bl;
		rect camera_rect = {camera_bl, camera_size};
		m3f screen_mat = m3f::translate(W/2, H/2) * m3f::scale(H, -H) * camera_mat.fast_inverse();

		v2f mouse_world;
		mouse_world.x = mouse_x;
		mouse_world.y = mouse_y;
		mouse_world = screen_mat.fast_inverse() * mouse_world;

		if ((mouse_state & SDL_BUTTON(SDL_BUTTON_RIGHT)) && !(prev_mouse_state & SDL_BUTTON(SDL_BUTTON_RIGHT)))
		{
			ball b(mouse_world, 0.0625f / zoom);
			b.is_moving = true;
			if (world_rect.contains(b.get_bounding_box()))
			{
				world.insert(std::move(b));
			}
		}
		if (mouse_state & SDL_BUTTON(SDL_BUTTON_LEFT))
		{
			v2f hs = {0.0625f/zoom, 0.0625f/zoom};
			rect mouse_rect = {mouse_world - hs, hs*2};
			v2f p1 = mouse_rect.pos;
			v2f p2 = p1 + mouse_rect.size;
			p1 = screen_mat * p1;
			p2 = screen_mat * p2;
			draw_rect(scr, p1.x, p1.y, p2.x, p2.y, {1.0f, 1.0f, 1.0f});
			auto to_delete = world.get_intersect(mouse_rect);
			for (auto o : to_delete)
			{
				if (ball_move == o) ball_moving = false;
				world.detach_obj(o);
				delete o;
			}
		}

		if ((mouse_state & SDL_BUTTON(SDL_BUTTON_MIDDLE)) && !(prev_mouse_state & SDL_BUTTON(SDL_BUTTON_MIDDLE)))
		{
			auto balls_move = world.get_contains(mouse_world);
			if (!balls_move.empty())
			{
				ball_move = balls_move.front();
				ball_moving = true;
			}
		}

		if (!(mouse_state & SDL_BUTTON(SDL_BUTTON_MIDDLE)) && (prev_mouse_state & SDL_BUTTON(SDL_BUTTON_MIDDLE)) && ball_moving)
		{
			ball_move->o.v = ball_move->o.p - mouse_world;
			ball_move->o.is_moving = true;
			ball_moving = false;
			ball_move = nullptr;
		}

		constexpr uint64_t tick_div = 1;
		for (uint64_t i = 0; i < tick_div; i++)
		{
			std::queue<obj_with_area<ball> *> to_update;
			std::stack<quad_tree_node<ball> *> s;
			s.push(world.root);
			while (!s.empty())
			{
				quad_tree_node<ball> *p = s.top();
				s.pop();
				for (uint64_t i = 0; i < 4; i++)
					if (p->n[i] != nullptr) s.push(p->n[i]);
				for (auto o : p->l)
				{
					if (o->o.is_moving)
						to_update.push(o);
				}
			}
			while (!to_update.empty())
			{
				auto o = to_update.front();
				ball &b = o->o;
				to_update.pop();

				b.p += b.v * dt / tick_div;
				b.v *= pow(0.75f, dt / tick_div);
				if (b.v.length2() < 0.001f)
				{
					b.v = {0.0f, 0.0f};
					b.is_moving = false;
				}

				// collision

				auto maybe_colliding = world.get_intersect(b.get_bounding_box());
				maybe_colliding.remove(o);
				for (auto c : maybe_colliding)
				{
					v2f rel = c->o.p - b.p;
					float d = b.r + c->o.r;
					float l = rel.length();
					if (l < d)
					{
						v2f axis = ~rel;
						float m1 = b.r * b.r;
						float m2 = c->o.r * c->o.r;
						float im = 1.0f / (m1 + m2);

						// static
						float dist = l - d;
						b.p += axis * dist * (m2 * im);
						c->o.p -= axis * dist * (m1 * im);

						// dynamic
						v2f n1 = axis * axis.dot(b.v);
						v2f n2 = axis * axis.dot(c->o.v);
						v2f t1 = b.v - n1;
						v2f t2 = c->o.v - n2;
						v2f nn1 = (n1 * (m1 - m2) + n2 * (m2 * 2.0f)) * im;
						v2f nn2 = (n1 * (m1 * 2.0f) + n2 * (m2 - m1)) * im;
						b.v = t1 + nn1;
						c->o.v = t2 + nn2;
						c->o.is_moving = true;
						//to_update.push(c);
					}
				}

				if (b.p.x < -world_size + b.r + 1.0f)
				{
					b.p.x += -world_size + b.r + 1.0f - b.p.x;
					b.v.x *= -1.0f;
				}
				if (b.p.x > world_size - b.r - 1.0f)
				{
					b.p.x += world_size - b.r - 1.0f - b.p.x;
					b.v.x *= -1.0f;
				}
				if (b.p.y < -world_size + b.r + 1.0f)
				{
					b.p.y += -world_size + b.r + 1.0f - b.p.y;
					b.v.y *= -1.0f;
				}
				if (b.p.y > world_size - b.r - 1.0f)
				{
					b.p.y += world_size - b.r - 1.0f - b.p.y;
					b.v.y *= -1.0f;
				}

				world.update_obj(o);
			}
		}

		if (keys[SDLK_g])
		{
			constexpr rgb col_lookup[6] = {
				{0.0f, 0.0f, 1.0f},
				{1.0f, 0.0f, 0.0f},
				{1.0f, 0.0f, 1.0f},
				{0.0f, 1.0f, 0.0f},
				{1.0f, 1.0f, 0.0f},
				{1.0f, 1.0f, 1.0f}
			};
			std::stack<quad_tree_node<ball> *> to_draw;
			to_draw.push(world.root);
			while (!to_draw.empty())
			{
				quad_tree_node<ball> *p = to_draw.top();
				to_draw.pop();
				if (camera_rect.intersects(p->node_rect))
				{
					for (uint64_t i = 0; i < 4; i++)
						if (p->n[i] != nullptr)
							to_draw.push(p->n[i]);
					v2f p1 = p->node_rect.pos;
					v2f p2 = p1 + p->node_rect.size;
					p1 = screen_mat * p1;
					p2 = screen_mat * p2;
					int d = 0;
					draw_rect(scr, p1.x+d, p1.y-d, p2.x-d, p2.y+d, col_lookup[p->len % 6]);
				}
			}
		}
		auto balls_draw = world.get_intersect(camera_rect);
		for (auto q : balls_draw)
		{
			v2f p = screen_mat * q->o.p;
			float r = q->o.r * zoom * H;
			if (ball_moving && q == ball_move)
				draw_circle(scr, p.x, p.y, r, {0.0f, 1.0f, 0.0f});
			else draw_circle(scr, p.x, p.y, r, {1.0f, 1.0f, 1.0f});
		}
		if (ball_moving)
		{
			v2f bm_scr = screen_mat * ball_move->o.p;
			draw_line(scr, mouse_x, mouse_y, bm_scr.x, bm_scr.y, {0.0f, 1.0f, 0.0f});
		}

		if (SDL_MUSTLOCK(surface)) SDL_LockSurface(surface);
		for (uint64_t i = 0; i < W*H; i++)
		{
			uint8_t r, g, b;
			scr->pixels[i].to_3u8(r, g, b);
			((uint32_t *)surface->pixels)[i] = SDL_MapRGB(surface->format, r, g, b);
		}
		if (SDL_MUSTLOCK(surface)) SDL_UnlockSurface(surface);
		SDL_Flip(surface);
		prev_mouse_state = mouse_state;
	}
	delete scr;
	SDL_Quit();
	return 0;
}