class grid {
public:

  grid(int nx, int ny, float sig);
  ~grid();
  void set(int index, float value);
  float get(int index);
  void g_set(int index, float value);
  float g_get(int index);
  void d_set(int index, float value);
  float d_get(int index);
  void smooth(int k);
  grid& operator= (grid & g);
  grid(const grid &obj);

//private:
  const int nx;
  const int ny;
  const float sigma;
  float * parameter;
  //float * grid_150;
  float * diff;
};

