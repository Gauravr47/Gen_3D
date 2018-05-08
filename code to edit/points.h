#ifndef points_hpp
#define points_hpp

/*Structure to save Integers pair in 2D, IP# denotes interger point where # is the dimension */
struct Ip2  
{ int x, z;
  Ip2 (int x = 0, int z = 0): x(x), z(z) {}
  inline  bool operator==(const Ip2& p) const { return x==p.x&&z==p.z; }
  inline  bool operator!=(const Ip2& p) const { return x!=p.x||z!=p.z; }
  inline const Ip2& operator+ (const Ip2& p) const { return Ip2(x + p.x, z + p.z); }
  
};
/*Structure to save floating pair in 2D, IP# denotes floating point where # is the dimension */
struct Dp2
{
  double x;
  double z;
  Dp2 (double x = 0, double z = 0) : x(x), z(z) {}
  Dp2 (Ip2 ip) : x((double)ip.x), z((double)ip.z) {}
  inline  bool operator==(const Dp2& p) const { return x==p.x&&z==p.z; }
  inline  bool operator!=(const Dp2& p) const { return x!=p.x||z!=p.z; }

};

struct  Ip3
{ int x, y, z;
  Ip3 (int x = 0, int y = 0, int z =0): x(x), y(y), z(z) {}
  inline  bool operator==(const Ip3& p) const { return x==p.x&&y==p.y&&z==p.z; }
  inline  bool operator!=(const Ip3& p) const { return x!=p.x||y!=p.y||z!=p.z; }

};

 struct Dp3
{
  double x;
  double y;
  double z;
  Dp3(double x = 0, double y = 0, double z =0) : x(x), y(y), z(z) {}
  Dp3 (Ip3 i) : x((double)i.x), y((double)i.y), z((double)i.z) {}
  inline  bool operator==(const Dp3& p) const { return x==p.x&&y==p.y&&z==p.z; }
  inline  bool operator!=(const Dp3& p) const { return x!=p.x||y!=p.y||z!=p.z; }

};

//Functions to convert between the types

   inline Dp3 int2ToDouble3 ( const Ip3& p, const double z){ return Dp3( (double)p.x, (double)p.y, (double)z);}
   inline Ip3 double2ToInt3 (const Dp3& p, const double  z){ return Ip3((int)p.x, (int)p.y, (int)z);}
   inline Ip2 add (Ip2& p, Ip2& q) { return Ip2(p.x+q.x, p.z+q.z);}
   inline Ip2 add (Ip2& p, int x, int z) { return Ip2(p.x+x, p.z+z);}
   inline Ip3 add (Ip3& p, Ip3& q) { return Ip3(p.x+q.x, p.y+q.y, p.z+q.z);}
   inline Ip3 add (Ip3& p, int x, int y, int z) { return Ip3(p.x+x, p.y+y, p.z+z);}
   inline Dp2 add (Dp2& p, Dp2& q) { return Dp2(p.x+q.x, p.z+q.z);}
   inline Dp2 add (Dp2& p, double x, double z) { return Dp2(p.x+x, p.z+z);}
   inline Dp3 add (Dp3& p, Dp3& q) { return Dp3(p.x+q.x, p.y+q.y, p.z+q.z);}
   inline Dp3 add (Dp3& p, double x, double y, double z) { return Dp3(p.x+x, p.y+y, p.z+z);}



 #endif
