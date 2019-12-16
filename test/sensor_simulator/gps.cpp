#include "gps.h"

namespace sensor_simulator
{
namespace sensor
{

Gps::Gps(std::shared_ptr<Ekf> ekf):Sensor(ekf)
{
}

Gps::~Gps()
{
}

void Gps::send(uint32_t time)
{
	_gps_data.time_usec = time;
	_ekf->setGpsData(time, _gps_data);
}

void Gps::setData(const gps_message& gps)
{
	_gps_data = gps;
}

void Gps::stepHeight(float hgt_change)
{
	_gps_data.alt += hgt_change * 1e3f;
}

void Gps::stepHorizontalPosition(Vector2f hpos_change)
{
	float hposN_curr;
	float hposE_curr;
	map_projection_global_project((float)_gps_data.lat, (float)_gps_data.lon, &hposN_curr, &hposE_curr);
	Vector2f hpos_new = Vector2f{hposN_curr, hposE_curr} + hpos_change;
	double lat_new;
	double lon_new;
	map_projection_global_reproject(hpos_new(0), hpos_new(1), &lat_new, &lon_new);
	_gps_data.lon = (uint32_t)lon_new;
	_gps_data.lat = (uint32_t)lat_new;
}

} // namespace sensor
} // namespace sensor_simulator
