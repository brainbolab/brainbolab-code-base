import sys
sys.path.insert(1, "../");
import camera_calibration_lib as ccl
import numpy as np


cam1_r_images_dir = "./camera_1/";
cam2_r_images_dir = "./camera_2/";

n_cb_rows = 3;
n_cb_cols = 6;
world_scale = 1.;

cam1_r = ccl.CameraCalibrator("1r");
cam1_r.load_images(cam1_r_images_dir);
cam1_r.view_images();
cam1_r.init_corners(n_cb_rows, n_cb_cols, world_scale);
cam1_r.find_corners_auto();
cam1_r.view_images_w_grid();
cam1_r.calibrate();

cam2_r = ccl.CameraCalibrator("2r");
cam2_r.load_images(cam2_r_images_dir);
cam2_r.view_images();
cam2_r.init_corners(n_cb_rows, n_cb_cols, world_scale);
cam2_r.find_corners_auto();
cam2_r.view_images_w_grid();
cam2_r.calibrate();

stereo12_r = ccl.StereoCalibrator(cam1_r, cam2_r);
stereo12_r.calibrate();
stereo12_r.compute_projection_matrices();

cam1_coords = cam1_r.get_image_coords("r_19.png");
cam2_coords = cam2_r.get_image_coords("r_19.png");
coords_3d = stereo12_r.triangulate(cam1_coords, cam2_coords);

# Print distances between chessboard corners to check how well calibration and triangulation work.
# Distances are in units of chessboard blocks.
print(
  "Distance between lower left and upper left corners of chessboard:",
  np.sqrt(((coords_3d[:,0] - coords_3d[:,1])**2).sum() ));

print(
  "Distance between upper left and upper right corners of chessboard:",
  np.sqrt(((coords_3d[:,1] - coords_3d[:,2])**2).sum() ));

print(
  "Distance between lower right and upper right corners of chessboard:",
  np.sqrt(((coords_3d[:,2] - coords_3d[:,3])**2).sum() ));

print(
  "Distance between lower left and lower right corners of chessboard:",
  np.sqrt(((coords_3d[:,0] - coords_3d[:,3])**2).sum() ));

print(
  "Distance between lower left and upper right corners of chessboard:",
  np.sqrt(((coords_3d[:,0] - coords_3d[:,2])**2).sum() ));

print(
  "Distance between upper left and lower right corners of chessboard:",
  np.sqrt(((coords_3d[:,1] - coords_3d[:,3])**2).sum() ));

cam1_r_intrinsic_mat = np.array([
  [1974.487773269557, 0.0, 779.5933149717231],
  [0.0, 1968.8720223601222, 616.713174788812],
  [0.0, 0.0, 1.0]],
  dtype=np.float64);

cam1_r_dist_coef = np.array([
  [0.2568755396038332, -7.215009846074793, 0.011745577464620556, 0.007967616557864175, 43.35462352950956]],
  dtype=np.float64);

cam2_r_intrinsic_mat = np.array([
  [2765.1665892850733, 0.0, 950.7200680374605],
  [0.0, 2733.990351416927, 499.85109889073374],
  [0.0, 0.0, 1.0]],
  dtype=np.float64);

cam2_r_dist_coef = np.array([
  [-0.4503370702129634, 15.938229248589565, -0.03603978325894311, 0.038540995013659964, -121.41849296787524]],
  dtype=np.float64);

stereo12_R = np.array([
  [0.9638682345378128, 0.26614756156705377, -0.01111314180919279],
  [-0.1947608119489957, 0.6756475348532966, -0.7110334976466519],
  [-0.18173126472329362, 0.6875070065959872, 0.703070311777917]],
  dtype=np.float64);

stereo12_T = np.array([
  [-0.40152334318322613],
  [12.747525918345012],
  [15.204533088333939]],
  dtype=np.float64);

cam1 = ccl.CameraCalibrator("1");
cam1.set_intrinsic_mat(cam1_r_intrinsic_mat);
cam1.set_dist_coef(cam1_r_dist_coef);

cam2 = ccl.CameraCalibrator("2");
cam2.set_intrinsic_mat(cam2_r_intrinsic_mat);
cam2.set_dist_coef(cam2_r_dist_coef);

assert (cam1.intrinsic_mat == cam1_r.intrinsic_mat).all();
assert (cam1.dist_coef == cam1_r.dist_coef).all();

assert (cam2.intrinsic_mat == cam2_r.intrinsic_mat).all();
assert (cam2.dist_coef == cam2_r.dist_coef).all();

stereo12 = ccl.StereoCalibrator(cam1, cam2);
stereo12.set_extrinsic_params(stereo12_R, stereo12_T);
stereo12.compute_projection_matrices();

assert (stereo12.R == stereo12_r.R).all();
assert (stereo12.T == stereo12_r.T).all();
assert (stereo12.proj1 == stereo12_r.proj1).all();
assert (stereo12.proj2 == stereo12_r.proj2).all();

coords_3d_man = stereo12.triangulate(cam1_coords, cam2_coords);

assert (coords_3d_man == coords_3d).all();

cam1_r.find_corners_man();
cam1_r.view_images_w_grid();
cam1_r.calibrate();

cam2_r.find_corners_man();
cam2_r.view_images_w_grid();
cam2_r.calibrate();

stereo12_r = ccl.StereoCalibrator(cam1_r, cam2_r);
stereo12_r.calibrate();
stereo12_r.compute_projection_matrices();

cam1_coords = cam1_r.get_image_coords("r_19.png");
cam2_coords = cam2_r.get_image_coords("r_19.png");
coords_3d = stereo12_r.triangulate(cam1_coords, cam2_coords);

# Print distances between chessboard corners to check how well calibration and triangulation work.
# Distances are in units of chessboard blocks.
print(
  "Distance between lower left and upper left corners of chessboard:",
  np.sqrt(((coords_3d[:,0] - coords_3d[:,1])**2).sum() ));

print(
  "Distance between upper left and upper right corners of chessboard:",
  np.sqrt(((coords_3d[:,1] - coords_3d[:,2])**2).sum() ));

print(
  "Distance between lower right and upper right corners of chessboard:",
  np.sqrt(((coords_3d[:,2] - coords_3d[:,3])**2).sum() ));

print(
  "Distance between lower left and lower right corners of chessboard:",
  np.sqrt(((coords_3d[:,0] - coords_3d[:,3])**2).sum() ));

print(
  "Distance between lower left and upper right corners of chessboard:",
  np.sqrt(((coords_3d[:,0] - coords_3d[:,2])**2).sum() ));

print(
  "Distance between upper left and lower right corners of chessboard:",
  np.sqrt(((coords_3d[:,1] - coords_3d[:,3])**2).sum() ));
