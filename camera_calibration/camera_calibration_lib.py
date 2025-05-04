import cv2
import os
import numpy as np
import pandas as pd


class CameraCalibrator:
  def __init__(self, cam_id):
    if not isinstance(cam_id, str):
      err_msg = "The provided camera ID must be a string variable."
      raise AssertionError(err_msg);

    self.cam_id = cam_id;
    self.images_dir = None;
    self.cam_info = None;
    self.img_height = None;
    self.img_width = None;
    self.n_cb_rows = None;
    self.n_cb_cols = None;
    self.world_scale = None;
    self.intrinsic_mat = None;
    self.dist_coef = None;


  def load_images(self, images_dir):
    if not os.path.isdir(images_dir):
      err_msg = "The provided directory:\n  " + images_dir + "\ndoes not exist.";
      raise AssertionError(err_msg);

    if not images_dir[-1] == "/":
      images_dir = images_dir + "/";

    cam_info = pd.DataFrame(sorted(os.listdir(images_dir) ), columns=["filename"]);
    cam_info["images"] = list(
      map(
        lambda fname_fx: cv2.imread(images_dir + fname_fx),
        cam_info["filename"]) );

    # Image dimensions. Images should be the same size across the provided directory.
    img_height, img_width = cam_info.loc[0,"images"].shape[0:2];

    img_dim_match_ixs = list(
      map(
        lambda img_ix: img_ix.shape[0:2] == (img_height, img_width),
        cam_info["images"]) );

    if not all(img_dim_match_ixs):
      err_msg = "Not all images have the same dimensions.";
      raise AssertionError(err_msg);

    self.images_dir = images_dir;
    self.cam_info = cam_info;
    self.img_height = img_height;
    self.img_width = img_width;


  def init_corners(self, n_cb_rows, n_cb_cols, world_scale=1.0):
    if not isinstance(self.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the CameraCalibrator object.";

      raise AssertionError(err_msg);

    expected_vars = {"filename", "images"};
    missing_vars = expected_vars.difference(set(self.cam_info.columns) );

    if len(missing_vars) > 0:
      err_msg = \
        "The 'cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars);

      raise AssertionError(err_msg);

    n_images = len(self.cam_info["images"]);

    # Coordinates of squares in the chessboard world space. 3D point in real world space.
    objp_temp = np.zeros((n_cb_rows * n_cb_cols, 3), np.float32);
    objp_temp[:,:2] = np.mgrid[0:n_cb_rows,0:n_cb_cols].T.reshape(-1, 2);
    objp_temp = world_scale * objp_temp;

    # Pixel coordinates of chessboards. 2D points in image plane.
    corners_temp = np.empty((n_cb_rows * n_cb_cols, 1, 2), dtype=np.float32);
    corners_temp[:] = np.nan;

    image_w_grid_temp = np.zeros((self.img_height, self.img_width, 3), dtype=np.uint8);

    self.cam_info["objpoints"] = [objp_temp.copy() for ix in range(0, n_images)];
    self.cam_info["corners"] = [corners_temp.copy() for ix in range(0, n_images)];
    self.cam_info["images_w_grid"] = [image_w_grid_temp.copy() for ix in range(0, n_images)];
    self.cam_info["found_corners"] = False;

    self.n_cb_rows = n_cb_rows;
    self.n_cb_cols = n_cb_cols;
    self.world_scale = world_scale;


  def find_corners_auto(self):
    if not isinstance(self.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the CameraCalibrator object.";

      raise AssertionError(err_msg);

    expected_vars = {"filename", "images", "objpoints", "corners", "images_w_grid", "found_corners"};
    missing_vars = expected_vars.difference(set(self.cam_info.columns) );

    if len(missing_vars) > 0:
      err_msg = \
        "The 'cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars);

      raise AssertionError(err_msg);

    if not isinstance(self.n_cb_rows, int) or isinstance(self.n_cb_rows, bool):
      err_msg = "The 'n_cb_rows' variable must be an integer.";
      raise AssertionError(err_msg);

    if not isinstance(self.n_cb_cols, int) or isinstance(self.n_cb_cols, bool):
      err_msg = "The 'n_cb_cols' variable must be an integer.";
      raise AssertionError(err_msg);

    if not isinstance(self.world_scale, float) or isinstance(self.world_scale, bool):
      err_msg = "The 'world_scale' variable must be a floating point number.";
      raise AssertionError(err_msg);

    # Criteria used by chessboard pattern detector. Change this if the code can't find the chessboard.
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.0001);

    n_images = len(self.cam_info["images"]);

    for ix in range(0, n_images):
      if not self.cam_info.loc[ix,"found_corners"]:
        img_ix = self.cam_info.loc[ix,"images"].copy();

        # Try to find the chessboard through simple grayscaling of the image.
        img_gray = cv2.cvtColor(img_ix, cv2.COLOR_BGR2GRAY);
        ret_ix, corners_ix = cv2.findChessboardCorners(img_gray, (self.n_cb_rows, self.n_cb_cols), None);

        if ret_ix:
          # Convolution size used to improve corner detection. Don't make this too large.
          conv_size = (11, 11);

          # Opencv can attempt to improve the chessboard coordinates.
          corners_ix = cv2.cornerSubPix(img_gray, corners_ix, conv_size, (-1, -1), criteria);
          cv2.drawChessboardCorners(img_ix, (self.n_cb_rows, self.n_cb_cols), corners_ix, ret_ix);
          self.cam_info.loc[ix,"images_w_grid"][:] = img_ix;
          self.cam_info.loc[ix,"corners"][:] = corners_ix;
          self.cam_info.loc[ix,"found_corners"] = ret_ix;

        else:
          warn_msg = "Cannot find chessboard corners in image {}.";
          print(warn_msg.format(self.cam_info.loc[ix,"filename"]) );


  def find_corners_man(self, refine_corners=True):
    if not isinstance(self.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the CameraCalibrator object.";

      raise AssertionError(err_msg);

    expected_vars = {"filename", "images", "objpoints", "corners", "images_w_grid", "found_corners"};
    missing_vars = expected_vars.difference(set(self.cam_info.columns) );

    if len(missing_vars) > 0:
      err_msg = \
        "The 'cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars);

      raise AssertionError(err_msg);

    # Criteria used by chessboard pattern detector. Change this if the code can't find the chessboard.
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.0001);

    n_images = len(self.cam_info["images"]);
    n_grid_points = self.n_cb_rows * self.n_cb_cols;

    def mouse_coords(event, x, y, flags, param):
      coords = param[0];
      img = param[1];
      img_title = param[2];

      if event == cv2.EVENT_LBUTTONDBLCLK:
        coords.append((x, y) );
        cv2.circle(img, (x, y), 5, (255, 0, 0), -1);
        cv2.imshow(img_title, img);

    for ix in range(0, n_images):
      if not self.cam_info.loc[ix,"found_corners"]:
        img_ix = self.cam_info.loc[ix,"images"].copy();
        img_title_ix = \
          "Cam ID: " + self.cam_id + " - " + \
          "Filename: " + self.cam_info.loc[ix,"filename"];

        coords_ix = [];

        cv2.imshow(img_title_ix, img_ix);
        cv2.setMouseCallback(img_title_ix, mouse_coords, (coords_ix, img_ix, img_title_ix) );
        cv2.waitKey(0);
        cv2.destroyAllWindows();

        if len(coords_ix) == n_grid_points:
          img_ix = self.cam_info.loc[ix,"images"].copy();
          corners_ix = np.array(coords_ix, dtype=np.float32).reshape((n_grid_points, 1, 2) );

          if refine_corners:
            img_gray = cv2.cvtColor(img_ix, cv2.COLOR_BGR2GRAY);

            # Convolution size used to improve corner detection. Don't make this too large.
            conv_size = (11, 11);

            # Opencv can attempt to improve the chessboard coordinates.
            corners_ix = cv2.cornerSubPix(img_gray, corners_ix, conv_size, (-1, -1), criteria);

          cv2.drawChessboardCorners(img_ix, (self.n_cb_rows, self.n_cb_cols), corners_ix, True);

          self.cam_info.loc[ix,"images_w_grid"][:] = img_ix;
          self.cam_info.loc[ix,"corners"][:] = corners_ix;
          self.cam_info.loc[ix,"found_corners"] = True;


  def calibrate(self):
    if not isinstance(self.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the CameraCalibrator object.";

      raise AssertionError(err_msg);

    expected_vars = {"corners", "found_corners", "objpoints"};
    missing_vars = expected_vars.difference(set(self.cam_info.columns) );

    if len(missing_vars) > 0:
      err_msg = \
        "The 'cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars);

      raise AssertionError(err_msg);

    rmse, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(
      list(self.cam_info.loc[self.cam_info["found_corners"],"objpoints"]),
      list(self.cam_info.loc[self.cam_info["found_corners"],"corners"]),
      (self.img_width, self.img_height),
      None,
      None);

    self.intrinsic_mat = mtx;
    self.dist_coef = dist;

    print("RMSE:", rmse);


  def undistort_points(self, coords):
    if (not isinstance(self.intrinsic_mat, np.ndarray) ) or (self.intrinsic_mat.shape != (3, 3) ):
      err_msg = "The 'intrinsic_mat' variable must be a 3x3 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(self.dist_coef, np.ndarray) ) or (self.dist_coef.shape != (1, 5) ):
      err_msg = "The 'dist_coef' variable must be a 1x5 Numpy array.";
      raise AssertionError(err_msg);

    n_points = coords.shape[0];

    coords_new = cv2.undistortPoints(
      coords,
      self.intrinsic_mat,
      self.dist_coef);

    coords_new = coords_new.reshape((n_points, 2) );

    return (
      np.concatenate([coords_new, np.ones([n_points, 1])], axis=1) @ \
      self.intrinsic_mat[0:2,:].T);


  def view_images(self):
    if not isinstance(self.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the CameraCalibrator object.";

      raise AssertionError(err_msg);

    expected_vars = {"filename", "images"};
    missing_vars = expected_vars.difference(set(self.cam_info.columns) );

    if len(missing_vars) > 0:
      err_msg = \
        "The 'cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars);

      raise AssertionError(err_msg);

    n_images = len(self.cam_info["images"]);

    for ix in range(0, n_images):
      img_ix = self.cam_info.loc[ix,"images"];
      img_title_ix = \
        "Cam ID: " + self.cam_id + " - " + \
        "Filename: " + self.cam_info.loc[ix,"filename"];

      cv2.imshow(img_title_ix, img_ix);
      cv2.waitKey(0);
      cv2.destroyAllWindows();


  def view_images_w_grid(self):
    if not isinstance(self.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the CameraCalibrator object.";

      raise AssertionError(err_msg);

    expected_vars = {"filename", "images_w_grid"};
    missing_vars = expected_vars.difference(set(self.cam_info.columns) );

    if len(missing_vars) > 0:
      err_msg = \
        "The 'cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars);

      raise AssertionError(err_msg);

    n_images = len(self.cam_info["images_w_grid"]);

    for ix in range(0, n_images):
      if self.cam_info.loc[ix,"found_corners"]:
        img_ix = self.cam_info.loc[ix,"images_w_grid"];
        img_title_ix = \
          "Cam ID: " + self.cam_id + " - " + \
          "Filename: " + self.cam_info.loc[ix,"filename"];

        cv2.imshow(img_title_ix, img_ix);
        cv2.waitKey(0);
        cv2.destroyAllWindows();


  def get_image_coords(self, img_name):
    if not isinstance(self.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the CameraCalibrator object.";

      raise AssertionError(err_msg);

    if not isinstance(img_name, str):
      err_msg = "The provided image name must be a string variable."
      raise AssertionError(err_msg);

    expected_vars = {"filename", "images"};
    missing_vars = expected_vars.difference(set(self.cam_info.columns) );

    if len(missing_vars) > 0:
      err_msg = \
        "The 'cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars);

      raise AssertionError(err_msg);

    ix_img = np.array(self.cam_info["filename"] == img_name);

    if not ix_img.any():
      err_msg = "Cannot find image with filename '" + img_name + "'";
      raise AssertionError(err_msg);

    ix_img = np.where(ix_img)[0][0];

    img = self.cam_info.loc[ix_img,"images"].copy();
    img_title = \
      "Cam ID: " + self.cam_id + " - " + \
      "Filename: " + self.cam_info.loc[ix_img,"filename"];

    coords = [];

    def mouse_coords(event, x, y, flags, param):
      coords = param[0];
      img = param[1];
      img_title = param[2];

      if event == cv2.EVENT_LBUTTONDBLCLK:
        coords.append((x, y) );
        cv2.circle(img, (x, y), 5, (255, 0, 0), -1);
        cv2.imshow(img_title, img);

    cv2.imshow(img_title, img);
    cv2.setMouseCallback(img_title, mouse_coords, (coords, img, img_title) );
    cv2.waitKey(0);
    cv2.destroyAllWindows();

    return np.array(coords, dtype=np.float64);


  def set_intrinsic_mat(self, intrinsic_mat):
    if (not isinstance(intrinsic_mat, np.ndarray) ) or (intrinsic_mat.shape != (3, 3) ):
      err_msg = "The provided 'intrinsic_mat' variable must be a 3x3 Numpy array.";
      raise AssertionError(err_msg);

    self.intrinsic_mat = intrinsic_mat.copy();


  def set_dist_coef(self, dist_coef):
    if (not isinstance(dist_coef, np.ndarray) ) or (dist_coef.shape != (1, 5) ):
      err_msg = "The provided 'dist_coef' variable must be a 1x5 Numpy array.";
      raise AssertionError(err_msg);

    self.dist_coef = dist_coef.copy();


class StereoCalibrator:
  def __init__(self, cam1, cam2):
    if not isinstance(cam1, CameraCalibrator):
      err_msg = "The provided 'cam1' variable must be a CameraCalibrator object."
      raise AssertionError(err_msg);

    if not isinstance(cam2, CameraCalibrator):
      err_msg = "The provided 'cam2' variable must be a CameraCalibrator object."
      raise AssertionError(err_msg);

    if (not isinstance(cam1.intrinsic_mat, np.ndarray) ) or (cam1.intrinsic_mat.shape != (3, 3) ):
      err_msg = "The 'cam1.intrinsic_mat' variable must be a 3x3 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(cam1.dist_coef, np.ndarray) ) or (cam1.dist_coef.shape != (1, 5) ):
      err_msg = "The 'cam1.dist_coef' variable must be a 1x5 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(cam2.intrinsic_mat, np.ndarray) ) or (cam2.intrinsic_mat.shape != (3, 3) ):
      err_msg = "The 'cam2.intrinsic_mat' variable must be a 3x3 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(cam2.dist_coef, np.ndarray) ) or (cam2.dist_coef.shape != (1, 5) ):
      err_msg = "The 'cam2.dist_coef' variable must be a 1x5 Numpy array.";
      raise AssertionError(err_msg);

    self.cam1 = cam1;
    self.cam2 = cam2;
    self.R = None;
    self.T = None;
    self.proj1 = None;
    self.proj2 = None;


  def calibrate(self):
    if not isinstance(self.cam1.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam1.cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the 'cam1' CameraCalibrator object.";

      raise AssertionError(err_msg);

    if not isinstance(self.cam2.cam_info, pd.DataFrame):
      err_msg = \
        "The 'cam2.cam_info' variable must be the Pandas DataFrame " + \
        "created by calling the 'load_images()' method on the 'cam2' CameraCalibrator object.";

      raise AssertionError(err_msg);

    expected_vars = {"corners", "filename", "found_corners", "objpoints"};
    missing_vars1 = expected_vars.difference(set(self.cam1.cam_info.columns) );

    if len(missing_vars1) > 0:
      err_msg = \
        "The 'cam1.cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars1);

      raise AssertionError(err_msg);

    missing_vars2 = expected_vars.difference(set(self.cam2.cam_info.columns) );

    if len(missing_vars2) > 0:
      err_msg = \
        "The 'cam2.cam_info' DataFrame is missing the following variables:\n  " + \
        "\n  ".join(missing_vars2);

      raise AssertionError(err_msg);

    if not isinstance(self.cam1.img_height, int) or isinstance(self.cam1.img_height, bool):
      err_msg = "The 'cam1.img_height' variable must be an integer.";
      raise AssertionError(err_msg);

    if not isinstance(self.cam1.img_width, int) or isinstance(self.cam1.img_width, bool):
      err_msg = "The 'cam1.img_width' variable must be an integer.";
      raise AssertionError(err_msg);

    if not isinstance(self.cam2.img_height, int) or isinstance(self.cam2.img_height, bool):
      err_msg = "The 'cam2.img_height' variable must be an integer.";
      raise AssertionError(err_msg);

    if not isinstance(self.cam2.img_width, int) or isinstance(self.cam2.img_width, bool):
      err_msg = "The 'cam2.img_width' variable must be an integer.";
      raise AssertionError(err_msg);

    if (self.cam1.img_height, self.cam1.img_width) != (self.cam2.img_height, self.cam2.img_width):
      err_msg = \
        "The image dimensions for 'cam1' do not match the image dimensions for 'cam2'.\n" + \
        "  cam1: " + str((self.cam1.img_height, self.cam1.img_width) ) + "\n" + \
        "  cam2: " + str((self.cam2.img_height, self.cam2.img_width) ) + "\n";

      raise AssertionError(err_msg);

    # Criteria used by chessboard pattern detector. Change this if the code can't find the chessboard.
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.0001);

    groups_to_keep = \
      self.cam1.cam_info[["filename", "found_corners"]].\
      merge(
        self.cam2.cam_info[["filename", "found_corners"]],
        how="inner",
        on="filename").\
      melt("filename")[["filename", "value"]].\
      groupby("filename").\
      all().\
      reset_index();

    groups_to_keep = pd.DataFrame(groups_to_keep.loc[groups_to_keep["value"],"filename"]);

    if len(groups_to_keep) > 0:
      cam_info1 = \
        self.cam1.cam_info[["filename", "objpoints", "corners"]].\
        merge(groups_to_keep, on="filename").\
        rename(columns={"corners": "corners1"});

      cam_info2 = \
        self.cam2.cam_info[["filename", "corners"]].\
        merge(groups_to_keep, on="filename").\
        rename(columns={"corners": "corners2"});

      stereo_calib_info = cam_info1.merge(cam_info2);

      rmse, CM1, dist1, CM2, dist2, R, T, E, F = cv2.stereoCalibrate(
        list(stereo_calib_info["objpoints"]),
        list(stereo_calib_info["corners1"]),
        list(stereo_calib_info["corners2"]),
        self.cam1.intrinsic_mat,
        self.cam1.dist_coef,
        self.cam2.intrinsic_mat,
        self.cam2.dist_coef,
        (self.cam1.img_width, self.cam1.img_height),
        criteria=criteria,
        flags=cv2.CALIB_FIX_INTRINSIC);

      self.R = R;
      self.T = T;

      print("RMSE:", rmse);

    else:
      print("No common image files with found chessboard coordinates between cameras.");


  def set_extrinsic_params(self, R, T):
    if (not isinstance(R, np.ndarray) ) or (R.shape != (3, 3) ):
      err_msg = "The provided extrinsic parameter rotation matrix 'R' variable must be a 3x3 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(T, np.ndarray) ) or (T.shape != (3, 1) ):
      err_msg = "The provided extrinsic parameter translation vector 'T' variable must be a 3x1 Numpy array.";
      raise AssertionError(err_msg);

    self.R = R.copy();
    self.T = T.copy();


  def compute_projection_matrices(self):
    if (not isinstance(self.R, np.ndarray) ) or (self.R.shape != (3, 3) ):
      err_msg = "The extrinsic parameter rotation matrix 'R' variable must be a 3x3 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(self.T, np.ndarray) ) or (self.T.shape != (3, 1) ):
      err_msg = "The extrinsic parameter translation vector 'T' variable must be a 3x1 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(self.cam1.intrinsic_mat, np.ndarray) ) or (self.cam1.intrinsic_mat.shape != (3, 3) ):
      err_msg = "The 'cam1.intrinsic_mat' variable must be a 3x3 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(self.cam2.intrinsic_mat, np.ndarray) ) or (self.cam2.intrinsic_mat.shape != (3, 3) ):
      err_msg = "The 'cam2.intrinsic_mat' variable must be a 3x3 Numpy array.";
      raise AssertionError(err_msg);

    # RT matrix for camera 1 is identity.
    RT1 = np.concatenate([np.eye(3), [[0],[0],[0]]], axis=-1);
    self.proj1 = self.cam1.intrinsic_mat @ RT1; # Projection matrix for camera 1.

    # RT matrix for camera 2 is the R and T obtained from stereo calibration.
    RT2 = np.concatenate([self.R, self.T], axis=-1);
    self.proj2 = self.cam2.intrinsic_mat @ RT2; # Projection matrix for camera 2.


  def triangulate(self, coords1, coords2):
    if (not isinstance(self.proj1, np.ndarray) ) or (self.proj1.shape != (3, 4) ):
      err_msg = "The 'proj1' variable must be a 3x4 Numpy array.";
      raise AssertionError(err_msg);

    if (not isinstance(self.proj2, np.ndarray) ) or (self.proj2.shape != (3, 4) ):
      err_msg = "The 'proj2' variable must be a 3x4 Numpy array.";
      raise AssertionError(err_msg);

    coords1_new = self.cam1.undistort_points(coords1);
    coords2_new = self.cam2.undistort_points(coords2);

    points_4d = cv2.triangulatePoints(self.proj1, self.proj2, coords1_new.T, coords2_new.T);

    return np.apply_along_axis(lambda col_ix: col_ix[0:3] / col_ix[3], 0, points_4d);
