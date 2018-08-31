#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// From miles per hour to meters per second
double mph2mps(double mph) { return mph / 2.23694; }
double mps2mph(double mps) { return mps * 2.23694; }

// Get average velocity within range
double get_avg_speed(vector<double> traffic)
{
    return accumulate(traffic.begin(), traffic.end(), 0.0) / traffic.size();
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s)
{
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("}");
    if (found_null != string::npos)
    {
        return "";
    }
    else if (b1 != string::npos && b2 != string::npos)
    {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{
    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for (int i = 0; i < maps_x.size(); i++)
    {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x, y, map_x, map_y);
        if (dist < closestLen)
        {
            closestLen = dist;
            closestWaypoint = i;
        }
    }

    return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double heading = atan2((map_y - y), (map_x - x));

    double angle = fabs(theta - heading);
    angle = min(2 * pi() - angle, angle);

    if (angle > pi() / 4)
    {
        closestWaypoint++;
        if (closestWaypoint == maps_x.size())
        {
            closestWaypoint = 0;
        }
    }

    return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

    int prev_wp;
    prev_wp = next_wp - 1;
    if (next_wp == 0)
    {
        prev_wp = maps_x.size() - 1;
    }

    double n_x = maps_x[next_wp] - maps_x[prev_wp];
    double n_y = maps_y[next_wp] - maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
    double proj_x = proj_norm * n_x;
    double proj_y = proj_norm * n_y;

    double frenet_d = distance(x_x, x_y, proj_x, proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000 - maps_x[prev_wp];
    double center_y = 2000 - maps_y[prev_wp];
    double centerToPos = distance(center_x, center_y, x_x, x_y);
    double centerToRef = distance(center_x, center_y, proj_x, proj_y);

    if (centerToPos <= centerToRef)
    {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for (int i = 0; i < prev_wp; i++)
    {
        frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
    }

    frenet_s += distance(0, 0, proj_x, proj_y);

    return {frenet_s, frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
    {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();

    double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s - maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
    double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

    double perp_heading = heading - pi() / 2;

    double x = seg_x + d * cos(perp_heading);
    double y = seg_y + d * sin(perp_heading);

    return {x, y};
}

int main()
{
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;

    ifstream in_map_(map_file_.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line))
    {
        istringstream iss(line);
        double x;
        double y;
        float s;
        float d_x;
        float d_y;
        iss >> x;
        iss >> y;
        iss >> s;
        iss >> d_x;
        iss >> d_y;
        map_waypoints_x.push_back(x);
        map_waypoints_y.push_back(y);
        map_waypoints_s.push_back(s);
        map_waypoints_dx.push_back(d_x);
        map_waypoints_dy.push_back(d_y);
    }

    // Set planner started flag
    bool _start = false;

    // Set speed limit
    float SPEED_LIMIT = 50; // MPH    

    h.onMessage([&_start, &SPEED_LIMIT, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                                                                                                                                    uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        //auto sdata = string(data).substr(0, length);
        //cout << sdata << endl;
        if (length && length > 2 && data[0] == '4' && data[1] == '2')
        {

            auto s = hasData(data);

            if (s != "")
            {
                auto j = json::parse(s);

                string event = j[0].get<string>();

                if (event == "telemetry")
                {
                    // j[1] is the data JSON object

                    // Main car's localization Data
                    double car_x = j[1]["x"];
                    double car_y = j[1]["y"];
                    double car_s = j[1]["s"];
                    double car_d = j[1]["d"];
                    double car_yaw = j[1]["yaw"];
                    double car_speed = j[1]["speed"];

                    // Previous path data given to the Planner
                    auto previous_path_x = j[1]["previous_path_x"];
                    auto previous_path_y = j[1]["previous_path_y"];
                    // Previous path's end s and d values
                    double end_path_s = j[1]["end_path_s"];
                    double end_path_d = j[1]["end_path_d"];

                    // Sensor Fusion Data, a list of all other cars on the same side of the road.
                    auto sensor_fusion = j[1]["sensor_fusion"];

                    // Previous path size
                    int prev_size = previous_path_x.size();                    

                    // Time for car to travel between each point
                    double _t = 0.02;

                    // Maximun for acceleration (or de-acceleration)
                    double MAX_ACC = 0.224;

                    // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

                    // Set reference speed
                    double ref_vel;
                    if (!_start)
                    {
                        ref_vel = car_speed;
                        _start = true;
                    }

                    if (prev_size > 0)
                    {
                        car_s = end_path_s;
                        car_d = end_path_d;
                    }

                    // Find current lane
                    int lane = car_d / 4;

                    // Keep track cars ahead, to the left or right
                    bool car_ahead = false;
                    bool car_left = false;
                    bool car_right = false;

                    // Flag if we want to change lane
                    bool lane_change = false;                    

                    // Loop through sensor fusion list
                    for (auto car : sensor_fusion)
                    {
                        double curr_d = car[6];
                        int car_lane = -1;

                        // Get car velocity
                        double curr_vx = car[3];
                        double curr_vy = car[4];
                        double curr_speed = sqrt(curr_vx * curr_vx + curr_vy * curr_vy);
                        double curr_s = car[5];
                        curr_s += (double)prev_size * curr_speed * _t;

                        // Determine which lane the current car is at
                        if (curr_d > lane * 4 && curr_d < (lane + 1) * 4)
                        {
                            // car in the same lane
                            car_lane = lane;                            
                        }
                        else if (curr_d > (lane - 1) * 4 && curr_d < lane * 4)
                        {
                            // car in the left lane
                            car_lane = lane - 1;                                
                        }
                        else if (curr_d > (lane + 1) * 4 && curr_d < (lane + 2) * 4)
                        {
                            // car in the right lane
                            car_lane = lane + 1;                            
                        }
                        else
                        {
                            // If not same or adjacent lane, skip this car
                            continue;
                        }

                        // Check if this car is within 30m ahead of ego car or 15m behind     
                        if (curr_s > car_s - 15 && curr_s < car_s + 30)
                        {
                            if (car_lane == lane && curr_s > car_s)
                            {
                                // Record car ahead
                                car_ahead = true;

                                if (mps2mph(curr_speed) < SPEED_LIMIT - 5)
                                {
                                    // If ahead car is slow, try to change lane
                                    lane_change = true;
                                }                                
                            }
                            else if (car_lane == lane - 1)
                            {
                                // Record car to the left
                                car_left = true;
                            }
                            else if (car_lane == lane + 1)
                            {
                                // Record car to the right
                                car_right = true;
                            }
                        }
                    }                    

                    // Adjust speed or change lane if needed
                    if (car_ahead) // If there is a car ahead
                    {
                        if (!car_left && lane > 0 && lane_change) // If no car to the left and there is a left lane
                        {
                            lane--;
                        }
                        else if (!car_right && lane < 2 && lane_change) // If no car to the right and there is a right lane
                        {
                            lane++;
                        }
                        else
                        {
                            // If not able to change lane, reduce speed
                            ref_vel -= MAX_ACC;
                        }
                    }
                    else if (ref_vel < SPEED_LIMIT - 0.5)
                    {
                        // If not too close to ahead vehicle, try to drive near speed limit
                        ref_vel += MAX_ACC;
                    }

                    // Anchor points for spline fit later
                    vector<double> ptsx;
                    vector<double> ptsy;

                    // Starting point, default is current car state
                    double ref_x = car_x;
                    double ref_y = car_y;
                    double ref_yaw = deg2rad(car_yaw);
                    double prev_ref_x;
                    double prev_ref_y;

                    // If previous path almost empty, use current car state
                    if (prev_size < 2)
                    {
                        // Find car position 1m before
                        prev_ref_x = ref_x - cos(ref_yaw);
                        prev_ref_y = ref_y - sin(ref_yaw);
                    }
                    // Use previous path's end point as starting
                    else
                    {
                        ref_x = previous_path_x[prev_size - 1];
                        ref_y = previous_path_y[prev_size - 1];

                        prev_ref_x = previous_path_x[prev_size - 2];
                        prev_ref_y = previous_path_y[prev_size - 2];

                        ref_yaw = atan2(ref_y - prev_ref_y, ref_x - prev_ref_x);
                    }

                    // Store these 2 points
                    ptsx.push_back(prev_ref_x);
                    ptsx.push_back(ref_x);
                    ptsy.push_back(prev_ref_y);
                    ptsy.push_back(ref_y);

                    // 30, 60, 90 m ahead of starting points in Frenet
                    vector<double> next_wp30 = getXY(car_s + 30, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    vector<double> next_wp60 = getXY(car_s + 60, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    vector<double> next_wp90 = getXY(car_s + 90, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

                    // Store these 3 points
                    ptsx.push_back(next_wp30[0]);
                    ptsx.push_back(next_wp60[0]);
                    ptsx.push_back(next_wp90[0]);

                    ptsy.push_back(next_wp30[1]);
                    ptsy.push_back(next_wp60[1]);
                    ptsy.push_back(next_wp90[1]);

                    // Transform to local (car) coordinates for easier calculation
                    for (int i = 0; i < ptsx.size(); i++)
                    {
                        // Shift car reference angle to 0 degree
                        double del_x = ptsx[i] - ref_x;
                        double del_y = ptsy[i] - ref_y;

                        ptsx[i] = del_x * cos(0 - ref_yaw) - del_y * sin(0 - ref_yaw);
                        ptsy[i] = del_x * sin(0 - ref_yaw) + del_y * cos(0 - ref_yaw);
                    }

                    // Initiate a spline and use anchor points to fit
                    tk::spline s;
                    s.set_points(ptsx, ptsy);

                    // x, y values sending to simulator
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;

                    // Start previous path points have not passed
                    for (int i = 0; i < previous_path_x.size(); i++)
                    {
                        next_x_vals.push_back(previous_path_x[i]);
                        next_y_vals.push_back(previous_path_y[i]);
                    }

                    // Calculate how to break spline points at desired speed
                    double target_x = 30.0;
                    double target_y = s(target_x);
                    double target_dist = sqrt(target_x * target_x + target_y * target_y);

                    double x_addon = 0.0;

                    // Fill up rest path
                    for (int i = 1; i <= 50 - previous_path_x.size(); i++)
                    {
                        double N = target_dist / (_t * mph2mps(ref_vel));
                        double x_point = x_addon + target_x / N;
                        double y_point = s(x_point);

                        x_addon = x_point;

                        // Shift back to original heading
                        double del_x = x_point;
                        double del_y = y_point;

                        x_point = del_x * cos(ref_yaw) - del_y * sin(ref_yaw);
                        y_point = del_x * sin(ref_yaw) + del_y * cos(ref_yaw);

                        x_point += ref_x;
                        y_point += ref_y;

                        next_x_vals.push_back(x_point);
                        next_y_vals.push_back(y_point);
                    }

                    json msgJson;
                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;

                    auto msg = "42[\"control\"," + msgJson.dump() + "]";

                    //this_thread::sleep_for(chrono::milliseconds(1000));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }
            }
            else
            {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1)
        {
            res->end(s.data(), s.length());
        }
        else
        {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port))
    {
        std::cout << "Listening to port " << port << std::endl;
    }
    else
    {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
