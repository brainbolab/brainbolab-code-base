import sys
sys.path.insert(1, "../");
import camera_calibration_lib as ccl
import numpy as np


sync1_images_dir = "./sync1/";
sync2_images_dir = "./sync2/";

n_cb_rows = 4;
n_cb_cols = 7;
world_scale = 1.;

cam1 = ccl.CameraCalibrator("1");
cam1.load_images(sync1_images_dir);
cam1.view_images();
cam1.init_corners(n_cb_rows, n_cb_cols, world_scale);
cam1.find_corners_auto();
cam1.view_images_w_grid();
cam1.calibrate();

cam2 = ccl.CameraCalibrator("2");
cam2.load_images(sync2_images_dir);
cam2.view_images();
cam2.init_corners(n_cb_rows, n_cb_cols, world_scale);
cam2.find_corners_auto();
cam2.view_images_w_grid();
cam2.calibrate();

stereo12 = ccl.StereoCalibrator(cam1, cam2);
stereo12.calibrate();
stereo12.compute_projection_matrices();

# These are the pixel coordinates of the chessboard corners from the 1.png files for each camera angle.
cam1_coords = np.array([[630, 434], [645, 97], [1053, 176], [1058,455]], dtype=np.float64);
cam2_coords = np.array([[234, 532], [246, 110], [926, 136], [905, 564]], dtype=np.float64);

# These are the 3D coordinates of the chessboard corners w.r.t the camera 1 angle.
coords_3d = stereo12.triangulate(cam1_coords, cam2_coords);

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
