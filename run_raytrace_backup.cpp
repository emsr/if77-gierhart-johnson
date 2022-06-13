
/*
g++ -g -DQUIET=0 -o run_raytrace_backup run_raytrace_backup.cpp
*/

#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * A point for atmospheric ray-tracing.
 */
struct RayTracePt
{
  /// Height of layer boundary [km].
  double height;
  /// Refractivity [dimensionless (_not_ in M-units)].
  double refractivity;
  /// Index of refraction [dimensionless].
  double index_of_refrac;
  /// radius of layer boundary [km].
  double radius;
  /// Logarithmic ratio.
  double a;
};

/**
 * 
 */
struct RayTrace
{
  /// Refractivity of the surface [M-units/km]
  double refractivity_surf;
  /// Elevation of transmitter above mean sea level [km].
  double elev_surf;
  /// First kilometer difference in refractivity above sea level
  double delta_refractivity;
  /// Exponential constant for refractivity decay with height [1/km].
  double ce;
  /// The radius of the surface at the transmitter.
  double radius_surf;
  /// Radius of the Earth [km].
  static constexpr double radius_earth = 6370.0;
  /// Number of height points in 'standard' atmosphere model.
  static constexpr int num_heights = 25;
  /// Array of ray-tracing points.
  RayTracePt pt[num_heights];
};

/// 'Standard' atmosphere heights [km].
constexpr double height[RayTrace::num_heights]
  {0.0, 0.01, 0.02, 0.05, 0.1,
   0.2, 0.305, 0.5, 0.7, 1.0,
   1.524, 2.0, 3.048, 5.0, 7.0,
   10.0, 20.0, 30.480, 50.0, 70.0,
   90.0, 110.0, 225.0, 350.0, 475.0};


/**
 * Fill the arrays required for ray-tracing given the surface elevation and refractivity.
 * On input the RayTrace object should have heights assigned. On output, te remainder of the
 * RayTracePt elements are assigned.
 *
 * @param[in,out] rt       Ray-tracing structure.
 * @param[in] refrac_surf  Refractivity at the surface [M-units]
 * @param[in] elev_surf    Elevation above sea-level of the surface [km]
 */
void
raytrace_init(RayTrace& rt, double refrac_surf, double elev_surf)
{
  rt.refractivity_surf = refrac_surf;
  rt.elev_surf = elev_surf;
  // Magic model for the first kilometer difference in refractivity above sea level [M-units/km]
  rt.delta_refractivity = -7.32 * std::exp(0.005577 * rt.refractivity_surf);
  rt.ce = std::log(rt.refractivity_surf / (rt.refractivity_surf + rt.delta_refractivity));
  rt.radius_surf = rt.radius_earth + rt.elev_surf;
  for (int i = 0; i < rt.num_heights; ++i)
  {
    rt.pt[i].height = height[i];
    rt.pt[i].refractivity = 1.0e-6 * rt.refractivity_surf * std::exp(-rt.ce * rt.pt[i].height);
    rt.pt[i].index_of_refrac = 1.0 + rt.pt[i].refractivity;
    rt.pt[i].radius = rt.radius_earth + rt.elev_surf + rt.pt[i].height;
  }
  rt.pt[0].a = 0.0;
  for (int i = 1; i < rt.num_heights; ++i)
  {
    const auto im1 = i - 1;
    //const auto a_numer = std::log(rt.pt[i].index_of_refrac) - std::log(rt.pt[im1].index_of_refrac); // could use log1p
    const auto a_numer = std::log1p(rt.pt[i].refractivity) - std::log1p(rt.pt[im1].refractivity);
    const auto a_denom = std::log(rt.pt[i].radius) - std::log(rt.pt[im1].radius);
    rt.pt[i].a = a_numer / a_denom;
  }
}

/**
 * Print ray-tracing array structure for warm fuzzy feeling.
 */
void
raytrace_write(const RayTrace& rt)
{
  if constexpr (!QUIET)
  {
    std::cout.precision(16);
    std::cout << std::showpoint;
    //std::cout << std::scientific;
    std::cout << "ENS = " << rt.refractivity_surf << '\n';
    std::cout << "DN = " << rt.delta_refractivity << '\n';
    std::cout << "HS = " << rt.elev_surf << '\n';
    std::cout << "az = " << rt.radius_earth << '\n';
    std::cout << "as = " << rt.radius_surf << '\n';
    for (int i = 0; i < rt.num_heights; ++i)
    {
      std::cout << ' ' << std::setw(23) << rt.pt[i].height
                << ' ' << std::setw(23) << rt.pt[i].refractivity
                << ' ' << std::setw(23) << rt.pt[i].index_of_refrac
                << ' ' << std::setw(23) << rt.pt[i].radius
                << ' ' << std::setw(23) << rt.pt[i].a
                << '\n';
    }
  }
}

/**
 * Given a ray-tracer object, the heights above the surface and a takeoff angle,
 * compute the ray elevation angle at the receive antenna and the great circle distance along the surface
 * to the receiver.
 *
 * @param[in] rt           Ray-tracing structure.
 * @param[in] height_tmit  Transmitter antenna height [km]
 * @param[in] height_recv  Receiver antenna height [km]
 * @param[in] theta_init   Input initial elevation or takeoff angle [rad].
 * @return The output ray parameters at the receiver.
 */
double
raytrace(const RayTrace& rt, double height_tmit, double height_recv, double theta_init, double& distance)
{
  auto theta_curr = theta_init;
  const auto radius_tmit = rt.radius_earth + rt.elev_surf + height_tmit;
  const auto refractivity_tmit = 1.0e-6 * rt.refractivity_surf * std::exp(-rt.ce * height_tmit);
  const auto index_of_refrac_tmit = 1.0 + refractivity_tmit;
  const auto radius_recv = rt.radius_earth + rt.elev_surf + height_recv;
  const auto refractivity_recv = 1.0e-6 * rt.refractivity_surf * std::exp(-rt.ce * height_recv);
  const auto index_of_refrac_recv = 1.0 + refractivity_recv;
  auto theta_all = 0.0;
  double theta[rt.num_heights];
  double slant_range = 0.0;

  bool theta_neg = true;

  if (theta_curr >= 0.0)
  {
    theta_neg = false;
  }
  else
  {
    double teg;
    if (rt.pt[0].radius == radius_tmit)
    {
      teg = 0.0;
    }
    else
    {
      const auto x = rt.pt[0].radius / (2.0 * radius_tmit);
      const auto z = (radius_tmit - rt.pt[0].radius) / rt.pt[0].radius;
      const auto w = (rt.pt[0].refractivity - refractivity_tmit) / index_of_refrac_tmit;
      teg = -2.0 * std::asin(std::sqrt(x * (z - w)));
    }
    if (theta_curr < teg)
      theta_curr = teg;
  }
  auto ate = std::abs(theta_curr);

  if (theta_curr >= 0.0)
    theta_neg = false;

  int nla = 0;
  if (theta_neg)
  {
    double ct;
    auto i_save = rt.num_heights - 1;
    for (int i = 1; i < rt.num_heights; ++i)
    {
      const auto y = 2.0 * sin(0.5 * ate) * sin(0.5 * ate);
      const auto z = (rt.pt[i].radius - radius_tmit) / radius_tmit;
      const auto w = (refractivity_tmit - rt.pt[i].refractivity) * std::cos(ate) / rt.pt[i].index_of_refrac;
      const auto x = y + z - w;
      if (x < 0.0)
        continue;
      ct = std::sqrt(0.5 * radius_tmit * x / rt.pt[i].radius);
      if (ct <= 1.0)
      {
        i_save = i;
        break;
      }
    }
    ct = 2.0 * std::asin(ct);
    theta_all = -2.0 * ct * rt.pt[i_save].a / (rt.pt[i_save].a + 1.0);
    theta[i_save] = ct;
    auto nk = i_save + 1;
    for (int i = nk; i < rt.num_heights; ++i)
    {
      auto radius_curr = rt.pt[i].radius;
      auto n_curr = rt.pt[i].index_of_refrac;
      if (radius_curr - radius_tmit > 0.0)
      {
        radius_curr = radius_tmit;
        n_curr = index_of_refrac_tmit;
      }
      const auto im1 = i - 1;
      auto ratio = rt.pt[im1].index_of_refrac * rt.pt[im1].radius / (n_curr * radius_curr);
      theta[i] = std::acos(ratio * std::cos(theta[im1]));
      auto x = -2.0 * rt.pt[i].a / (rt.pt[i].a + 1.0);
      theta_all += x * (theta[i] - theta[im1]);
      nla = i;
      if (radius_curr == radius_tmit)
        break;
    }
  }
  else
  {
    nla = rt.num_heights;
    for (int nl = 1; nl < rt.num_heights; ++nl)
    {
      if (radius_tmit <= rt.pt[nl].radius)
      {
        nla = nl;
        break;
      }
    }
  }

  if (nla < 1)
    nla = 1;
  theta[nla - 1] = ate;
  double theta_final = 0.0;
  bool recv_above_layers = true;
  for (int i = nla; i < rt.num_heights; ++i)
  {
    const auto im1 = i - 1;
    auto radius_curr = rt.pt[i].radius;
    auto n_curr = rt.pt[i].index_of_refrac;
    auto refrac_curr = rt.pt[i].refractivity;
    if (radius_curr - radius_recv > 0.0)
    {
      n_curr = index_of_refrac_recv;
      radius_curr = radius_recv;
      refrac_curr = refractivity_recv;
    }
    auto x = radius_tmit / (2.0 * radius_curr);
    const auto y = 2.0 * std::sin(0.5 * ate) * std::sin(0.5 * ate);
    const auto z = (radius_curr - radius_tmit) / radius_tmit;
    const auto w = (refractivity_tmit - refrac_curr) * std::cos(ate) / n_curr;
    theta[i] = 2.0 * std::asin(std::sqrt(x * (y + z - w)));
    x = -rt.pt[i].a / (rt.pt[i].a + 1.0);
    theta_all += x * (theta[i] - theta[im1]);
    theta_final = theta[i];
    if (rt.pt[i].radius > radius_recv)
      recv_above_layers = false;
  }
  if (recv_above_layers)
  {
    const auto i_last = rt.num_heights - 1;
    const auto x = rt.pt[i_last].index_of_refrac * rt.pt[i_last].radius / radius_recv;
    theta_final = std::acos(x * std::cos(theta_final));
  }
  const auto ca = theta_final - theta_curr + theta_all;
  distance = rt.radius_surf * ca; // We need to return this too.
  const auto ct = std::cos(theta_all);
  const auto st = std::sin(theta_all);
  const auto tnt = std::tan(theta_final);
  const auto y = index_of_refrac_recv / index_of_refrac_tmit;
  auto x = (ct - st * tnt - y) / (y * std::tan(theta_curr) - st - ct * tnt);
  x = std::atan(x);
  auto cx = theta_curr - x; // Unused?

  return theta_final;
}

/**
 * Ray trace from start height to specified end height and distance.
 * Iteratively solve d(\theta) - d_F = 0 for theta;
 *
 * @param[in] rt           Ray-tracing structure.
 * @param[in] height_tmit  Transmitter antenna height [km]
 * @param[in] height_recv  Receiver antenna height [km]
 * @param[in] distance     Great-circle sea-level distance from transmitter to receiver [km].
 * @param[in] theta_init   Input initial elevation or takeoff angle [rad].
 * @param[in] tol          The distance tolerance required for success [km].
 * @param[out] theta_final Output takeoff angle such that the ray arrives at the given receiver height and distance.
 * @return The output ray parameters at the receiver.
 */
double
raytrace_point(const RayTrace& rt, double height_tmit, double height_recv, double distance, double theta_init, double tol,
               double& distance_final)
{
  double dist_0, dist_1;
  auto theta_0 = theta_init;
  auto theta_1 = theta_0 + 0.005;
  raytrace(rt, height_tmit, height_recv, theta_0, dist_0);
  raytrace(rt, height_tmit, height_recv, theta_1, dist_1);
  for (int i = 0; i < 10; ++i)
  {
    const auto jacobian_inv = (theta_1 - theta_0) / (dist_1 - dist_0);
    const auto theta_2 = theta_1 - jacobian_inv * (dist_1 - distance);
    double dist_2;
    raytrace(rt, height_tmit, height_recv, theta_2, dist_2);
    theta_0 = theta_1;
    dist_0 = dist_1;
    theta_1 = theta_2;
    dist_1 = dist_2;
    if constexpr (!QUIET)
      std::cout << ' ' << i+1 << ' ' << theta_2 << ' ' << dist_2 << '\n';
    if (std::abs(dist_1 - distance) < tol)
      break;
  }

  distance_final = dist_1;
  return theta_1;
}

/**
 * Return the approximate takeoff angle assuming an effective earth radius
 * based on the delta-refractivity at the surface.
 *
 * @param[in] rt           Ray-tracing structure.
 * @param[in] height_tmit  Transmit antenna height [km]
 * @param[in] height_recv  Receive antenna height [km]
 * @param[in] distance     Great-circle sea-level distance from transmit to receive [km].
 * @return takeoff angle of ray in radians above horizontal.
 */
double
takeoff_angle(const RayTrace& rt, double height_tmit, double height_recv, double distance)
{
  auto k = 1.0 / (1.0 + rt.radius_earth * 1.0e-6 * rt.delta_refractivity);
  auto radius_eff = k * rt.radius_surf;
  auto correction = 1.0 / (2.0 * radius_eff);
  auto diff_height = height_recv - correction * distance * distance - height_tmit;
  //auto slant_range = std::sqrt(distance * distance + diff_height * diff_height);
  return std::atan2(diff_height, distance);
}

/*
 * Main test driver.
 */
int
main()
{
  const auto zero = 0.0;
  auto hrp = 0.0; // ground plane elev (km presumably)
  auto h1 = 4.8768 / 1000; // km
  auto h2 = 0.3048 * 30000 / 1000; // km
  auto hp2 = h2 - hrp;
  auto hp1 = h1 - hrp;

  RayTrace rt;

  if constexpr (!QUIET)
    std::cout << '\n';
  raytrace_init(rt, 329.0, hrp);
  raytrace_write(rt);
  double distance;
  auto ry = raytrace(rt, hp1, hp2, -1.56, distance);
  auto ds0 = distance;
  if constexpr (!QUIET)
  {
    std::cout << "ry = " << ry << '\n';
    std::cout << "qqd = " << distance << '\n';
  }

  if constexpr (!QUIET)
    std::cout << '\n';
  raytrace_init(rt, 301.0, hrp);
  raytrace_write(rt);
  auto ry2 = raytrace(rt, zero, hp2, zero, distance);
  auto horizon_recv = distance;
  if constexpr (!QUIET)
  {
    std::cout << "ry2 = " << ry2 << '\n';
    std::cout << "qqd = " << distance << '\n';
  }

  if constexpr (!QUIET)
    std::cout << '\n';
  auto ry1 = raytrace(rt, zero, hp1, zero, distance);
  auto horizon_tmit = distance;
  if constexpr (!QUIET)
  {
    std::cout << "ry1 = " << ry1 << '\n';
    std::cout << "qqd = " << distance << '\n';
  }

  // Made up....
  if constexpr (!QUIET)
    std::cout << '\n';
  auto ray = raytrace(rt, hp1, hp2, 0.05, distance);
  auto dist = distance;
  if constexpr (!QUIET)
  {
    std::cout << "ray = " << ray << '\n';
    std::cout << "qqd = " << distance << '\n';
  }

  // Made up point downward....
  if constexpr (!QUIET)
    std::cout << '\n';
  auto ray2 = raytrace(rt, hp1, hp2, -0.05, distance);
  auto dist2 = distance;
  if constexpr (!QUIET)
  {
    std::cout << "ray2 = " << ray2 << '\n';
    std::cout << "qqd2 = " << distance << '\n';
  }

  // Now try to trace to a point!
  if constexpr (!QUIET)
    std::cout << '\n';
  auto distance_f = 25.0;
  auto theta_0 = takeoff_angle(rt, hp1, hp2, distance_f);
  auto theta_f = raytrace_point(rt, hp1, hp2, distance_f, theta_0, 0.001, distance);
  if constexpr (!QUIET)
  {
    std::cout << "distance_f = " << distance_f << '\n';
    std::cout << "theta_0 = " << theta_0 << '\n';
    std::cout << "theta_f = " << theta_f << '\n';
    std::cout << "distance = " << distance << '\n';
  }

  if constexpr (!QUIET)
    std::cout << '\n';
  distance_f = 100.0;
  theta_0 = takeoff_angle(rt, hp1, hp2, distance_f);
  theta_f = raytrace_point(rt, hp1, hp2, distance_f, theta_0, 0.001, distance);
  if constexpr (!QUIET)
  {
    std::cout << "distance_f = " << distance_f << '\n';
    std::cout << "theta_0 = " << theta_0 << '\n';
    std::cout << "theta_f = " << theta_f << '\n';
    std::cout << "distance = " << distance << '\n';
  }
}

