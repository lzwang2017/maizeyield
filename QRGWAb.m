function bt = QRGWAb(x, y, W, QP_A, QP_b, btL, btU, x0, QP_Aeq, QP_beq)
  p = size(x, 2);
  QP_c = -x'*(W.^2.*y);
  QP_Q = x'*(repmat(W.^2,1,p).*x);

  QP_c = round(QP_c*1e6)/1e6;
  QP_Q = round(QP_Q*1e6)/1e6;
  QP_A = round(QP_A*1e6)/1e6;
  QP_b = round(QP_b*1e6)/1e6;

  dd = 1e-6*min([1; QP_Q(QP_Q(:)>0)]);
  mineg = min(eig(QP_Q));
  while mineg < 1e-6
    dd = dd + max(abs(mineg), 1e-6);
    QP_Q = QP_Q + dd*eye(p);
    QP_Q = round(QP_Q*1e6)/1e6;
    mineg = min(eig(QP_Q));
  end

  options.Display = 'on';
  if ~exist('x0', 'var')
    x0 = [];
  end
  if ~exist('QP_Aeq', 'var')
    QP_Aeq = zeros(0, p);
  end
  if ~exist('QP_beq', 'var')
    QP_beq = [];
  end

  if exist ('OCTAVE_VERSION', 'builtin')
    options.MaxIter = 2e3;
    [bt, obj, info, lambda] = qp(x0, QP_Q, QP_c, QP_Aeq, QP_beq, btL, btU, -inf(size(QP_b)), QP_A, QP_b, options);
  else
    bt = gquadprog(QP_Q, QP_c, QP_A, QP_b, QP_Aeq, QP_beq, btL, btU, x0, options);
end
%   if isempty(bt)
%       bt = gquadprog(QP_Q, QP_c, QP_A, QP_b, zeros(0, p), [], btL, btU, x0, options);
%   end
end
