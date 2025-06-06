program vicsek_simulation
  implicit none

  ! --- パラメータ定義 ---
  integer, parameter :: n_particles = 100     ! 粒子数 (N)
  real(8), parameter :: box_lx = 10.0d0       ! シミュレーションボックスのX方向の長さ
  real(8), parameter :: box_ly = 10.0d0       ! シミュレーションボックスのY方向の長さ
  real(8), parameter :: r0 = 1.0d0            ! 相互作用半径 (R0)
  real(8), parameter :: v0 = 0.1d0            ! 粒子の速さ
  real(8), parameter :: eta = 0.4d0           ! ノイズ強度、パラメータの変更はここをいじる
  real(8), parameter :: dt = 1.0d0            ! 時間ステップ (Δt)
  integer, parameter :: total_steps = 4000    ! 総シミュレーションステップ数
  integer, parameter :: save_interval = 10    ! データ保存間隔
  real(8), parameter :: PI = acos(-1.0d0)     ! 円周率

  ! --- 配列宣言 ---
  real(8), dimension(n_particles) :: pos_x, pos_y         ! 時刻 t+Δt の粒子座標
  real(8), dimension(n_particles) :: pos_x_p, pos_y_p     ! 時刻 t の粒子座標 (previous)
  real(8), dimension(n_particles) :: theta              ! 粒子の角度
  real(8), dimension(n_particles) :: vel_x, vel_y         ! 粒子の速度ベクトル成分

  ! --- その他変数 ---
  integer :: ts, i
  character(len=20) :: filename_format
  character(len=30) :: output_filename

  ! ディレクトリtrajが存在することを確認してください
  ! mkdir traj コマンドで事前に作成してください。
  filename_format = "(A,I4.4,A)" ! "traj/pos", step_num, ".d"

  ! --- 粒子の初期化 ---
  call initialize_particles(n_particles, pos_x_p, pos_y_p, theta, box_lx, box_ly, PI)

  ! 最初の座標を現座標にもコピー
  pos_x = pos_x_p
  pos_y = pos_y_p

  ! --- メインタイムステップループ (スライド12ページ参照) ---
  do ts = 1, total_steps
    ! 1. 角度の更新 (スライド11ページ参照)
    !    時刻tの座標(pos_x_p, pos_y_p)と角度(theta)を使用
    call update_angles_vicsek(n_particles, theta, pos_x_p, pos_y_p, box_lx, box_ly, r0, eta, PI)

    ! 2. 位置の更新 (スライド12ページ参照)
    !    更新された角度(theta)と時刻tの座標(pos_x_p, pos_y_p)から時刻t+Δtの座標(pos_x, pos_y)を計算
    call update_positions_vicsek(n_particles, pos_x, pos_y, pos_x_p, pos_y_p, &
                                 theta, vel_x, vel_y, v0, dt, box_lx, box_ly)

    ! 3. 次のステップのために座標を保存 (スライド12ページ参照)
    pos_x_p = pos_x
    pos_y_p = pos_y

    ! 4. ファイルへの保存 (スライド12ページ参照)
    if (mod(ts, save_interval) == 0) then
      write(output_filename, filename_format) 'traj2/pos', ts / save_interval, '.d' !ここの出力フォルダ名を変える
      call save_data_vicsek(output_filename, n_particles, pos_x, pos_y, vel_x, vel_y)
      print *, "Saved data at step: ", ts, " to ", trim(output_filename)
    end if

    if (mod(ts,100) == 0) then
        print *, "Completed step: ", ts
    end if

  end do

  print *, "Simulation finished."

contains

  ! サブルーチン: initialize_particles
  ! 目的: 粒子の初期位置と角度をランダムに設定する
  subroutine initialize_particles(n, p_x, p_y, th, lx, ly, pi_val)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: p_x, p_y, th
    real(8), intent(in) :: lx, ly, pi_val
    integer :: i_particle
    real(8), dimension(n) :: rand_nums_pos_x, rand_nums_pos_y, rand_nums_theta

    call random_seed() ! 乱数シードを初期化

    ! 位置をランダムに設定 (0 から lx, 0 から ly の範囲)
    call random_number(rand_nums_pos_x)
    p_x = rand_nums_pos_x * lx
    call random_number(rand_nums_pos_y)
    p_y = rand_nums_pos_y * ly

    ! 角度をランダムに設定 (-PI から PI の範囲)
    call random_number(rand_nums_theta)
    th = (rand_nums_theta * 2.0d0 - 1.0d0) * pi_val
  end subroutine initialize_particles

  ! サブルーチン: update_angles_vicsek (スライド11ページ参照)
  ! 目的: Vicsekモデルに基づいて粒子の角度を更新する
  subroutine update_angles_vicsek(n, th, p_x_current, p_y_current, lx, ly, r_interaction, noise_eta, pi_val)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(inout) :: th
    real(8), dimension(n), intent(in) :: p_x_current, p_y_current
    real(8), intent(in) :: lx, ly, r_interaction, noise_eta, pi_val

    real(8), dimension(n) :: sum_s_x, sum_s_y       ! 各粒子の近傍ベクトルの合計成分
    real(8), dimension(n) :: current_s_x, current_s_y ! 各粒子の現時刻のベクトル成分
    real(8) :: dx, dy, dist_sq
    real(8) :: r0_sq
    real(8) :: random_val
    integer :: i_particle, j_particle

    r0_sq = r_interaction**2

    ! 1. 現時刻の各粒子の速度ベクトル成分を計算 (cos(theta), sin(theta))
    do i_particle = 1, n
      current_s_x(i_particle) = cos(th(i_particle))
      current_s_y(i_particle) = sin(th(i_particle))
    end do

    ! 2. 各粒子の合計ベクトルを自身のベクトルで初期化
    sum_s_x = current_s_x
    sum_s_y = current_s_y

    ! 3. 近傍粒子の速度ベクトルを加算 (スライド11ページの二重ループと条件分岐)
    do i_particle = 1, n
      do j_particle = 1, n
        if (i_particle == j_particle) cycle ! 自分自身とのペアはスキップ (既に初期化で考慮済み)

        ! 粒子iと粒子jの距離ベクトルを計算
        dx = p_x_current(i_particle) - p_x_current(j_particle)
        dy = p_y_current(i_particle) - p_y_current(j_particle)

        ! 周期境界条件を適用 (スライド9ページ参照、dnintを使用)
        dx = dx - dnint(dx / lx) * lx
        dy = dy - dnint(dy / ly) * ly

        dist_sq = dx**2 + dy**2

        ! 距離が相互作用半径内か判定
        if (dist_sq < r0_sq) then
          sum_s_x(i_particle) = sum_s_x(i_particle) + current_s_x(j_particle)
          sum_s_y(i_particle) = sum_s_y(i_particle) + current_s_y(j_particle)
        end if
      end do
    end do
    ! スライド11の do i=1,n-1 と do j=i+1,n の形にする場合は以下のようになるが、
    ! 相互作用は対称なので、上記のように全ペアを調べて自身のベクトルを初期値とする方が考えやすい。
    ! もしスライドのループ構造に忠実にする場合は、
    ! sum_s_x = current_s_x (初期化)
    ! sum_s_y = current_s_y
    ! do i_particle = 1, n - 1
    !   do j_particle = i_particle + 1, n
    !     ... (dx, dy, dist_sq の計算) ...
    !     if (dist_sq < r0_sq) then
    !       sum_s_x(i_particle) = sum_s_x(i_particle) + current_s_x(j_particle)
    !       sum_s_y(i_particle) = sum_s_y(i_particle) + current_s_y(j_particle)
    !       sum_s_x(j_particle) = sum_s_x(j_particle) + current_s_x(i_particle) ! 対称性
    !       sum_s_y(j_particle) = sum_s_y(j_particle) + current_s_y(i_particle)
    !     end if
    !   end do
    ! end do
    ! 上記は同等のはず。現在の実装は、各粒子iについて他の全ての粒子jを調べて加算する形。

    ! 4. 新しい角度を計算 (スライド11ページ参照: atan2とノイズ)
    do i_particle = 1, n
      call random_number(random_val) ! 0から1までの一様乱数
      th(i_particle) = atan2(sum_s_y(i_particle), sum_s_x(i_particle)) + &
                       noise_eta * pi_val * (2.0d0 * random_val - 1.0d0)
      ! 角度を -PI から PI の範囲に正規化 (atan2は通常この範囲だが念のため)
      if (th(i_particle) > pi_val) then
        th(i_particle) = th(i_particle) - 2.0d0 * pi_val
      else if (th(i_particle) < -pi_val) then
        th(i_particle) = th(i_particle) + 2.0d0 * pi_val
      end if
    end do
  end subroutine update_angles_vicsek

  ! サブルーチン: update_positions_vicsek (スライド12ページ参照)
  ! 目的: 粒子の位置を更新し、周期境界条件を適用する
  subroutine update_positions_vicsek(n, p_x_new, p_y_new, p_x_old, p_y_old, th, v_x, v_y, speed, time_step, lx, ly)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: p_x_new, p_y_new  ! 更新後の座標
    real(8), dimension(n), intent(in) :: p_x_old, p_y_old    ! 更新前の座標
    real(8), dimension(n), intent(in) :: th                  ! 粒子の角度
    real(8), dimension(n), intent(out) :: v_x, v_y           ! 速度ベクトル成分 (保存用)
    real(8), intent(in) :: speed, time_step, lx, ly
    integer :: i_particle

    do i_particle = 1, n
      ! 速度ベクトル成分の計算 (スライド8ページ参照)
      v_x(i_particle) = speed * cos(th(i_particle))
      v_y(i_particle) = speed * sin(th(i_particle))

      ! 位置の更新 (r_new = r_old + v * dt)
      p_x_new(i_particle) = p_x_old(i_particle) + v_x(i_particle) * time_step
      p_y_new(i_particle) = p_y_old(i_particle) + v_y(i_particle) * time_step

      ! 周期境界条件の適用 (スライド5ページ参照、floorを使用)
      ! pos(i) = pos(i) - lx * floor(pos(i) / lx)
      p_x_new(i_particle) = p_x_new(i_particle) - lx * floor(p_x_new(i_particle) / lx)
      p_y_new(i_particle) = p_y_new(i_particle) - ly * floor(p_y_new(i_particle) / ly)

      ! 別の周期境界条件の適用方法 (0未満ならlxを足し、lx以上ならlxを引く)
      ! if (p_x_new(i_particle) < 0.0d0) p_x_new(i_particle) = p_x_new(i_particle) + lx
      ! if (p_x_new(i_particle) >= lx) p_x_new(i_particle) = p_x_new(i_particle) - lx
      ! if (p_y_new(i_particle) < 0.0d0) p_y_new(i_particle) = p_y_new(i_particle) + ly
      ! if (p_y_new(i_particle) >= ly) p_y_new(i_particle) = p_y_new(i_particle) - ly
    end do
  end subroutine update_positions_vicsek

  ! サブルーチン: save_data_vicsek (スライド12ページ参照)
  ! 目的: 粒子の座標と速度ベクトル成分をファイルに保存する
  subroutine save_data_vicsek(f_name, n, p_x, p_y, v_comp_x, v_comp_y)
    implicit none
    character(len=*), intent(in) :: f_name
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: p_x, p_y, v_comp_x, v_comp_y
    integer :: i_particle
    integer, parameter :: file_unit = 23 ! ファイルユニット番号 (任意)

    open(unit=file_unit, file=trim(f_name), status='replace', action='write')
    do i_particle = 1, n
      write(file_unit, *) p_x(i_particle), p_y(i_particle), v_comp_x(i_particle), v_comp_y(i_particle)
    end do
    close(file_unit)
  end subroutine save_data_vicsek

end program vicsek_simulation
